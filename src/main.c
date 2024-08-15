#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "graph.h"
#include "local_search.h"
#include "mwis.h"

static inline int compare(const void *a, const void *b)
{
    return *(int *)a - *(int *)b;
}

static inline int compare_r(const void *a, const void *b, void *c)
{
    return ((int *)c)[*(int *)a] - ((int *)c)[*(int *)b];
}

graph reduction(graph g, long long *offset)
{
    *offset = 0;
    int *mask = malloc(sizeof(int) * g.N);
    int *rm = malloc(sizeof(int) * g.N);

    for (int i = 0; i < g.N; i++)
        mask[i] = 1;

    int imp = 1;
    while (imp)
    {
        imp = 0;
        for (int u = 0; u < g.N; u++)
        {
            if (!mask[u])
                continue;
            int nw = 0;
            for (int i = g.V[u]; i < g.V[u + 1]; i++)
            {
                int v = g.E[i];
                if (!mask[v])
                    continue;
                nw += g.W[v];
            }

            if (nw <= g.W[u])
            {
                imp = 1;
                *offset += g.W[u];
                mask[u] = 0;
                for (int i = g.V[u]; i < g.V[u + 1]; i++)
                    mask[g.E[i]] = 0;
            }
        }
    }

    graph sg = graph_subgraph(g, mask, rm);
    free(mask);
    free(rm);
    return sg;
}

void explore(graph g, local_search *ls, int k, int verbose)
{
    int c = 0, t = 0;
    double t0 = omp_get_wtime();

    long long best = ls->c;

    int inc = 0, exc = 0;

    int ref = ls->lc;

    while (c++ < k)
    {

        int _inc = 0, _exc = 0;

        ls->lc = ref;
        int p = ls->lc;

        int u = rand() % g.N;
        while (ls->T[u])
            u = rand() % g.N;

        if (ls->IS[u])
            _exc++;
        else
            _inc++;

        if (ls->IS[u])
            local_search_remove_vertex_tabu(g, ls, u);
        else
            local_search_add_vertex_tabu(g, ls, u);

        int r = rand() % 4;
        for (int i = 0; i < r && ls->qc > 0; i++)
        {
            int v = ls->Q[rand() % ls->qc];
            int q = 0;
            while (q++ < 100 && ls->T[v])
                v = ls->Q[rand() % ls->qc];

            if (q == 100)
                continue;

            if (ls->IS[v])
                _exc++;
            else
                _inc++;

            if (ls->IS[v])
                local_search_remove_vertex(g, ls, v);
            else
                local_search_add_vertex(g, ls, v);
        }

        local_search_greedy(g, ls);

        if (ls->c > best)
        {
            t++;
            best = ls->c;
            inc += _inc;
            exc += _exc;
            if (verbose)
            {
                printf("\r%lld %.2lf %d %d %d %d    ", ls->c, omp_get_wtime() - t0, c, t, inc, exc);
                fflush(stdout);
            }
        }
        else if (ls->c < best)
        {
            local_search_unwind(g, ls, p);
        }
        local_search_unwind_tabu(g, ls, p);
    }
    if (verbose)
        printf("\n");
}

int *run_ls(graph g, int *IS, int rp, int k, long long offset)
{
    local_search *ls = local_search_init(g);

    for (int i = 0; i < g.N; i++)
        if (IS[i])
            local_search_add_vertex(g, ls, i);

    local_search_greedy(g, ls);

    for (int i = 0; i < g.N; i++)
        ls->P[i] = 0;

    explore(g, ls, k, 1);

    int *P = malloc(sizeof(int) * g.N);
    for (int i = 0; i < g.N; i++)
        P[i] = i;

    local_search *lsr = local_search_init(g);

    qsort_r(P, g.N, sizeof(int), compare_r, ls->P);

    for (int i = 0; i < g.N / 2; i++)
    {
        if (ls->IS[P[i]])
            local_search_add_vertex_tabu(g, lsr, P[i]);
    }

    int *S = malloc(sizeof(int) * g.N);
    for (int i = 0; i < g.N; i++)
        S[i] = ls->IS[i];

    qsort(ls->P, g.N, sizeof(int), compare);

    for (int i = 0; i < 30; i++)
        printf("%d ", ls->P[i]);
    printf("\n");

    for (int i = 0; i < 30; i++)
        printf("%d ", ls->P[g.N - 1 - i]);
    printf("\n");

    local_search_free(ls);

    return S;
}

int main(int argc, char **argv)
{
    FILE *f = fopen(argv[1], "r");
    graph g = graph_parse(f);
    fclose(f);

    int *IS = malloc(sizeof(int) * g.N);
    for (int i = 0; i < g.N; i++)
        IS[i] = 0;

    if (argc > 2)
    {
        f = fopen(argv[2], "r");
        int u;
        while (fscanf(f, "%d\n", &u) == 1)
            IS[u - 1] = 1;
        fclose(f);
    }

    long long is = mwis_validate(g, IS);
    int valid = graph_validate(g.N, g.V, g.E);
    printf("Graph valid %d |V|=%d, |E|=%d IS=%lld\n", valid, g.N, g.V[g.N] / 2, is);

    int *S = run_ls(g, IS, 8, 4000000, 0);

    long long val = mwis_validate(g, S);
    printf("%s %lld\n", argv[1], val);

    free(IS);
    free(S);
    graph_free(g);

    return 0;

    // int *C = malloc(sizeof(int) * g.N);
    // for (int i = 0; i < g.N; i++)
    //     C[i] = 0;

    // int it = 200;
    // for (int i = 0; i < it; i++)
    // {
    //     int *S = run_ls(g, 64, 100000, 0);
    //     int in = 0, out = 0;
    //     for (int j = 0; j < g.N; j++)
    //     {
    //         C[j] += S[j];
    //         if (C[j] == 0)
    //             out++;
    //         else if (C[j] == i + 1)
    //             in++;
    //     }
    //     printf("%d %d / %d\n", in, out, g.N);
    //     free(S);
    // }

    // int *mask = malloc(sizeof(int) * g.N);
    // int offset = 0;
    // for (int i = 0; i < g.N; i++)
    //     mask[i] = 1;

    // for (int u = 0; u < g.N; u++)
    // {
    //     if (C[u] < it)
    //         continue;
    //     offset += g.W[u];
    //     mask[u] = 0;
    //     for (int i = g.V[u]; i < g.V[u + 1]; i++)
    //     {
    //         int v = g.E[i];
    //         mask[v] = 0;
    //     }
    // }

    // graph sg = graph_subgraph(g, mask, C);

    // for (int i = 0; i < 10; i++)
    // {
    //     int *S = run_ls(sg, 8, 10000000, offset);
    //     free(S);
    // }

    // graph_free(g);

    // return 0;
}