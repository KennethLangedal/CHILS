#include "graph.h"

#include <stdlib.h>
#include <sys/mman.h>

static inline void parse_id(char *data, size_t *p, int *v)
{
    while (data[*p] < '0' || data[*p] > '9')
        (*p)++;

    *v = 0;
    while (data[*p] >= '0' && data[*p] <= '9')
        *v = (*v) * 10 + data[(*p)++] - '0';
}

graph graph_parse(FILE *f)
{
    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *data = mmap(0, size, PROT_READ, MAP_PRIVATE, fileno_unlocked(f), 0);
    size_t p = 0;

    int N, M, t;
    parse_id(data, &p, &N);
    parse_id(data, &p, &M);
    parse_id(data, &p, &t);

    int *V = malloc(sizeof(int) * (N + 1));
    int *E = malloc(sizeof(int) * (M * 2));

    int *W = malloc(sizeof(int) * N);

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
            parse_id(data, &p, E + ei);
            E[ei++]--;
        }
        p++;
    }
    V[N] = ei;

    munmap(data, size);

    return (graph){.N = N, .V = V, .E = E, .W = W};
}

void graph_free(graph g)
{
    free(g.V);
    free(g.E);
    free(g.W);
}

static inline int compare(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

int graph_validate(int N, const int *V, const int *E)
{
    int M = 0;
    for (int u = 0; u < N; u++)
    {
        if (V[u + 1] - V[u] < 0)
            return 0;

        M += V[u + 1] - V[u];

        for (int i = V[u]; i < V[u + 1]; i++)
        {
            if (i < 0 || i >= V[N])
                return 0;

            int v = E[i];
            if (v < 0 || v >= N || v == u || (i > V[u] && v <= E[i - 1]))
                return 0;

            if (bsearch(&u, E + V[v], V[v + 1] - V[v], sizeof(int), compare) == NULL)
                return 0;
        }
    }

    if (M != V[N])
        return 0;

    return 1;
}

graph graph_subgraph(graph g, int *mask, int *rm)
{
    int *fm = malloc(sizeof(int) * g.N);
    int N = 0, M = 0;
    for (int u = 0; u < g.N; u++)
    {
        if (!mask[u])
            continue;

        fm[u] = N;
        rm[N] = u;
        N++;

        for (int i = g.V[u]; i < g.V[u + 1]; i++)
            if (mask[g.E[i]])
                M++;
    }

    graph sg = {.N = N};
    sg.V = malloc(sizeof(int) * (g.N + 1));
    sg.E = malloc(sizeof(int) * M);
    sg.W = malloc(sizeof(int) * g.N);

    M = 0;
    for (int u = 0; u < g.N; u++)
    {
        if (!mask[u])
            continue;

        sg.W[fm[u]] = g.W[u];
        sg.V[fm[u]] = M;

        for (int i = g.V[u]; i < g.V[u + 1]; i++)
        {
            int v = g.E[i];
            if (!mask[v])
                continue;

            sg.E[M] = fm[v];
            M++;
        }
    }
    sg.V[sg.N] = M;
    free(fm);
    return sg;
}