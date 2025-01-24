#include "graph.h"

#include <omp.h>
#include <stdlib.h>
#include <sys/mman.h>

static inline void parse_id(char *data, size_t *p, long long *v)
{
    while (data[*p] < '0' || data[*p] > '9')
        (*p)++;

    *v = 0;
    while (data[*p] >= '0' && data[*p] <= '9')
        *v = (*v) * 10 + data[(*p)++] - '0';
}

graph *graph_parse(FILE *f)
{
    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *data = mmap(0, size, PROT_READ, MAP_PRIVATE, fileno_unlocked(f), 0);
    size_t p = 0;

    long long N, M, t;
    parse_id(data, &p, &N);
    parse_id(data, &p, &M);
    parse_id(data, &p, &t);

    int *V = malloc(sizeof(int) * (N + 1));
    int *E = malloc(sizeof(int) * (M * 2));

    long long *W = malloc(sizeof(long long) * N);

    int ei = 0;
    for (int u = 0; u < N; u++)
    {
        parse_id(data, &p, W + u);
        V[u] = ei;
        while (ei < M * 2)
        {
            while (data[p] == ' ')
                p++;
            if (data[p] == '\n')
                break;

            long long e;
            parse_id(data, &p, &e);
            E[ei++] = e - 1;
            ;
        }
        p++;
    }
    V[N] = ei;

    munmap(data, size);

    graph *g = malloc(sizeof(graph));
    *g = (graph){.N = N, .V = V, .E = E, .W = W};

    return g;
}

// store graph in metis format
void graph_store(FILE *f, graph *g)
{
    fprintf(f, "%d %d 10\n", g->N, g->V[g->N] / 2);
    for (int u = 0; u < g->N; u++)
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
    for (int u = 0; u < g->N; u++)
    {
        if (g->V[u + 1] - g->V[u] < 0)
            return 0;

        M += g->V[u + 1] - g->V[u];

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            if (i < 0 || i >= g->V[g->N])
                return 0;

            int v = g->E[i];
            if (v < 0 || v >= g->N || v == u || (i > g->V[u] && v <= g->E[i - 1]))
                return 0;

            // if (bsearch(&u, g->E + g->V[v], g->V[v + 1] - g->V[v], sizeof(int), compare) == NULL)
            //     return 0;
        }
    }

    if (M != g->V[g->N])
        return 0;

    return 1;
}

graph *graph_subgraph(graph *g, int *mask, int *reverse_map)
{
    int *forward_map = malloc(sizeof(int) * g->N);
    int N = 0, M = 0;
    for (int u = 0; u < g->N; u++)
    {
        if (!mask[u])
            continue;

        forward_map[u] = N;
        reverse_map[N] = u;
        N++;

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
            if (mask[g->E[i]])
                M++;
    }

    graph *sg = malloc(sizeof(graph));
    *sg = (graph){.N = N};

    sg->V = malloc(sizeof(int) * (N + 1));
    sg->E = malloc(sizeof(int) * M);
    sg->W = malloc(sizeof(long long) * N);

    M = 0;
    for (int u = 0; u < g->N; u++)
    {
        if (!mask[u])
            continue;

        sg->W[forward_map[u]] = g->W[u];
        sg->V[forward_map[u]] = M;

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (!mask[v])
                continue;

            sg->E[M] = forward_map[v];
            M++;
        }
    }
    sg->V[sg->N] = M;

    free(forward_map);

    return sg;
}

void graph_subgraph_par(graph *g, graph *sg, int *mask, int *reverse_map, int *forward_map, int *s1, int *s2)
{
    int nt = omp_get_num_threads();
    int tid = omp_get_thread_num();
    int N = 0, M = 0;
#pragma omp for nowait
    for (int u = 0; u < g->N; u++)
    {
        if (!mask[u])
            continue;

        N++;

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
            if (mask[g->E[i]])
                M++;
    }

    s1[tid] = N;
    s2[tid] = M;

#pragma omp barrier

    int No = 0, Mo = 0;
    for (int i = 0; i < tid; i++)
    {
        No += s1[i];
        Mo += s2[i];
    }

    N = 0;
#pragma omp for
    for (int u = 0; u < g->N; u++)
    {
        forward_map[u] = N + No;
        reverse_map[N + No] = u;
        N++;

        sg->W[forward_map[u]] = g->W[u];
        sg->V[forward_map[u]] = M + Mo;
    }

    if (tid == nt)
    {
        sg->N = No + N;
        sg->V[sg->N] = M + Mo;
    }

    M = 0;
#pragma omp for
    for (int u = 0; u < g->N; u++)
    {
        if (!mask[u])
            continue;

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (!mask[v])
                continue;

            sg->E[M] = forward_map[v];
            M++;
        }
    }
}