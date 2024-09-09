#include "clique_graph.h"

#include <stdlib.h>
#include <sys/mman.h>

static inline int compare(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

static inline void parse_id(char *data, size_t *p, long long *v)
{
    while (data[*p] < '0' || data[*p] > '9')
        (*p)++;

    *v = 0;
    while (data[*p] >= '0' && data[*p] <= '9')
        *v = (*v) * 10 + data[(*p)++] - '0';
}

clique_graph *clique_graph_parse_json(FILE *f)
{
    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *data = mmap(0, size, PROT_READ, MAP_PRIVATE, fileno_unlocked(f), 0);
    size_t p = 0;

    while (data[p] != '[')
        p++;

    size_t _p = p;
    int N = 1;
    while (data[_p] != ']')
        if (data[_p++] == ',')
            N++;

    long long *W = malloc(sizeof(long long) * N);

    for (int i = 0; i < N; i++)
        parse_id(data, &p, W + i);

    while (data[p] != '[')
        p++;
    p++;

    _p = p;
    int NC = 0, M = 0;
    while (data[_p] != '}')
    {
        if (data[_p++] != '[')
            continue;
        NC++, M++;
        while (data[_p] != ']')
            if (data[_p++] == ',')
                M++;
    }

    int *C = malloc(sizeof(int) * (NC + 1));
    int *EC = malloc(sizeof(int) * M);

    int *D = malloc(sizeof(int) * N);
    for (int i = 0; i < N; i++)
        D[i] = 0;

    NC = 0, M = 0;
    while (data[p] != '}')
    {
        if (data[p++] != '[')
            continue;
        C[NC++] = M;
        long long c;
        do
        {
            parse_id(data, &p, &c);
            EC[M++] = c;
            D[c]++;
        } while (data[p] == ',');
    }
    C[NC] = M;

    for (int i = 0; i < NC; i++)
        qsort(EC + C[i], C[i + 1] - C[i], sizeof(int), compare);

    int *V = malloc(sizeof(int) * (N + 1));
    V[0] = 0;
    for (int i = 1; i <= N; i++)
        V[i] = V[i - 1] + D[i - 1];

    for (int i = 0; i < N; i++)
        D[i] = 0;

    int *EV = malloc(sizeof(int) * V[N]);

    for (int c = 0; c < NC; c++)
        for (int i = C[c]; i < C[c + 1]; i++)
            EV[V[EC[i]] + D[EC[i]]++] = c;

    munmap(data, size);

    free(D);

    clique_graph *cg = malloc(sizeof(clique_graph));
    *cg = (clique_graph){.N = N, .NC = NC, .V = V, .EV = EV, .C = C, .EC = EC, .W = W};

    return cg;
}

void clique_graph_free(clique_graph *g)
{
    free(g->V);
    free(g->EV);
    free(g->C);
    free(g->EC);
    free(g->W);
    free(g);
}

int clique_graph_validate(clique_graph *g)
{
    for (int u = 0; u < g->N; u++)
    {
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int c = g->EV[i];
            if (c < 0 || c >= g->NC || (i > g->V[u] && c <= g->EV[i - 1]))
                return 0;

            if (bsearch(&u, g->EC + g->C[c], g->C[c + 1] - g->C[c], sizeof(int), compare) == NULL)
                return 0;
        }
    }

    for (int c = 0; c < g->NC; c++)
    {
        for (int i = g->C[c]; i < g->C[c + 1]; i++)
        {
            int u = g->EC[i];
            if (u < 0 || u >= g->N || (i > g->C[c] && u <= g->EC[i - 1]))
                return 0;

            if (bsearch(&c, g->EV + g->V[u], g->V[u + 1] - g->V[u], sizeof(int), compare) == NULL)
                return 0;
        }
    }

    return 1;
}