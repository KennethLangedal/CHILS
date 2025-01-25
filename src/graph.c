#include "graph.h"

#include <omp.h>
#include <stdlib.h>
#include <sys/mman.h>

static inline void parse_id(char *Data, size_t *p, long long *v)
{
    while (Data[*p] < '0' || Data[*p] > '9')
        (*p)++;

    *v = 0;
    while (Data[*p] >= '0' && Data[*p] <= '9')
        *v = (*v) * 10 + Data[(*p)++] - '0';
}

graph *graph_parse(FILE *f)
{
    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *Data = mmap(0, size, PROT_READ, MAP_PRIVATE, fileno_unlocked(f), 0);
    size_t p = 0;

    long long n, m, t;
    parse_id(Data, &p, &n);
    parse_id(Data, &p, &m);
    parse_id(Data, &p, &t);

    int *V = malloc(sizeof(int) * (n + 1));
    int *E = malloc(sizeof(int) * (m * 2));

    long long *W = malloc(sizeof(long long) * n);

    int ei = 0;
    for (int u = 0; u < n; u++)
    {
        parse_id(Data, &p, W + u);
        V[u] = ei;
        while (ei < m * 2)
        {
            while (Data[p] == ' ')
                p++;
            if (Data[p] == '\n')
                break;

            long long e;
            parse_id(Data, &p, &e);
            E[ei++] = e - 1;
            ;
        }
        p++;
    }
    V[n] = ei;

    munmap(Data, size);

    graph *g = malloc(sizeof(graph));
    *g = (graph){.n = n, .V = V, .E = E, .W = W};

    return g;
}

// store graph in metis format
void graph_store(FILE *f, graph *g)
{
    fprintf(f, "%d %d 10\n", g->n, g->V[g->n] / 2);
    for (int u = 0; u < g->n; u++)
    {
        fprintf(f, "%lld", g->W[u]);
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
            fprintf(f, " %d", g->E[i] + 1);
        fprintf(f, "\n");
    }
}

void graph_free(graph *g)
{
    if (g == NULL)
        return;

    free(g->V);
    free(g->E);
    free(g->W);

    free(g);
}

static inline int compare(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

int graph_validate(graph *g)
{
    int M = 0;
    for (int u = 0; u < g->n; u++)
    {
        if (g->V[u + 1] - g->V[u] < 0)
            return 0;

        M += g->V[u + 1] - g->V[u];

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            if (i < 0 || i >= g->V[g->n])
                return 0;

            int v = g->E[i];
            if (v < 0 || v >= g->n || v == u || (i > g->V[u] && v <= g->E[i - 1]))
                return 0;

            // if (bsearch(&u, g->E + g->V[v], g->V[v + 1] - g->V[v], sizeof(int), compare) == NULL)
            //     return 0;
        }
    }

    if (M != g->V[g->n])
        return 0;

    return 1;
}

graph *graph_subgraph(graph *g, int *Mask, int *RM)
{
    int *FM = malloc(sizeof(int) * g->n);
    int n = 0, m = 0;
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        FM[u] = n;
        RM[n] = u;
        n++;

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
            if (Mask[g->E[i]])
                m++;
    }

    graph *sg = malloc(sizeof(graph));
    *sg = (graph){.n = n};

    sg->V = malloc(sizeof(int) * (n + 1));
    sg->E = malloc(sizeof(int) * m);
    sg->W = malloc(sizeof(long long) * n);

    m = 0;
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        sg->W[FM[u]] = g->W[u];
        sg->V[FM[u]] = m;

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (!Mask[v])
                continue;

            sg->E[m] = FM[v];
            m++;
        }
    }
    sg->V[sg->n] = m;

    free(FM);

    return sg;
}

void graph_subgraph_par(graph *g, graph *sg, int *Mask, int *RM, int *FM, int *S1, int *S2)
{
    int nt = omp_get_num_threads();
    int tid = omp_get_thread_num();
    int n = 0, m = 0;
#pragma omp for nowait
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        n++;

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
            if (Mask[g->E[i]])
                m++;
    }

    S1[tid] = n;
    S2[tid] = m;

#pragma omp barrier

    int n_o = 0, m_o = 0;
    for (int i = 0; i < tid; i++)
    {
        n_o += S1[i];
        m_o += S2[i];
    }

    n = n_o;
#pragma omp for
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        FM[u] = n;
        RM[n] = u;
        sg->W[n] = g->W[u];
        n++;
    }

    m = m_o;
#pragma omp for nowait
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        sg->V[FM[u]] = m;

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (!Mask[v])
                continue;

            sg->E[m] = FM[v];
            m++;
        }
    }

    if (tid == nt - 1)
    {
        sg->n = n;
        sg->V[sg->n] = m;
    }
#pragma omp barrier
}