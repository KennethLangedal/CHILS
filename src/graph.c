#include "graph.h"

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

void graph_free(graph *g)
{
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

            if (bsearch(&u, g->E + g->V[v], g->V[v + 1] - g->V[v], sizeof(int), compare) == NULL)
                return 0;
        }
    }

    if (M != g->V[g->N])
        return 0;

    return 1;
}

graph *graph_subgraph(graph *g, int *mask, int *reverse_mapping)
{
    int *forward_mapping = malloc(sizeof(int) * g->N);
    int N = 0, M = 0;
    for (int u = 0; u < g->N; u++)
    {
        if (!mask[u])
            continue;

        forward_mapping[u] = N;
        reverse_mapping[N] = u;
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

        sg->W[forward_mapping[u]] = g->W[u];
        sg->V[forward_mapping[u]] = M;

        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (!mask[v])
                continue;

            sg->E[M] = forward_mapping[v];
            M++;
        }
    }
    sg->V[sg->N] = M;

    free(forward_mapping);

    return sg;
}