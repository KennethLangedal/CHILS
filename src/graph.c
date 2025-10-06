#include "graph.h"

#include <omp.h>
#include <stdlib.h>
#include <limits.h>

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

    char *Data = malloc(size);
    size_t red = fread(Data, 1, size, f);
    size_t p = 0;

    long long n, m, t;
    parse_id(Data, &p, &n);
    parse_id(Data, &p, &m);
    parse_id(Data, &p, &t);

    if (n >= INT_MAX)
    {
        fprintf(stderr, "Number of vertices must be less than %d, got %lld\n", INT_MAX, n);
        exit(1);
    }

    long long *V = malloc(sizeof(long long) * (n + 1));
    int *E = malloc(sizeof(int) * (m * 2));

    long long *W = malloc(sizeof(long long) * n);

    long long ei = 0;
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

            if (e > n)
            {
                fprintf(stderr, "Edge endpoint out of bounds, {%lld, %lld}\n", u + 1ll, e);
                exit(1);
            }

            E[ei++] = e - 1;
        }
        p++;
    }
    V[n] = ei;

    free(Data);

    graph *g = malloc(sizeof(graph));
    *g = (graph){.n = n, .m = ei, .V = V, .E = E, .W = W};

    return g;
}

// store graph in metis format
void graph_store(FILE *f, graph *g)
{
    fprintf(f, "%d %lld 10\n", g->n, g->m / 2);
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

int graph_validate(graph *g)
{
    int M = 0;
    for (int u = 0; u < g->n; u++)
    {
        int d_u = g->V[u + 1] - g->V[u];
        if (d_u < 0)
            return 0;

        M += d_u;

        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
        {
            if (i < 0 || i >= g->m)
                return 0;

            int v = g->E[i];
            if (v < 0 || v >= g->n || v == u || (i > g->V[u] && v <= g->E[i - 1]))
                return 0;
        }
    }

    if (M != g->V[g->n] || M != g->m)
        return 0;

    return 1;
}

graph *graph_subgraph(graph *g, int *Mask, int *RM)
{
    int *FM = malloc(sizeof(int) * g->n);
    long long n = 0, m = 0;
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        FM[u] = n;
        RM[n] = u;
        n++;

        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
            if (Mask[g->E[i]])
                m++;
    }

    graph *sg = malloc(sizeof(graph));
    *sg = (graph){.n = n};

    sg->V = malloc(sizeof(long long) * (n + 1));
    sg->E = malloc(sizeof(int) * m);
    sg->W = malloc(sizeof(long long) * n);

    m = 0;
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        sg->W[FM[u]] = g->W[u];
        sg->V[FM[u]] = m;

        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
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

void graph_subgraph_par(graph *g, graph *sg, int *Mask, int *RM, int *FM, long long *S1, long long *S2)
{
    int nt = omp_get_num_threads();
    int tid = omp_get_thread_num();
    long long n = 0, m = 0;

#pragma omp for nowait
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        n++;

        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
            if (Mask[g->E[i]])
                m++;
    }

    S1[tid] = n;
    S2[tid] = m;

#pragma omp barrier

    long long n_offset = 0;
    for (int i = 0; i < tid; i++)
        n_offset += S1[i];

    n = n_offset;
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

    long long m_offset = 0;
    for (int i = 0; i < tid; i++)
        m_offset += S2[i];

    m = m_offset;
#pragma omp for nowait
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        sg->V[FM[u]] = m;

        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
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
        sg->m = m;
        sg->V[sg->n] = m;
    }
#pragma omp barrier
}