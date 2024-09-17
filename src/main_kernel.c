#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "graph.h"
#include "local_search.h"
#include "pils.h"
#include "reductions.h"

long long mwis_validate(graph *g, int *independent_set)
{
    long long cost = 0;
    for (int u = 0; u < g->N; u++)
    {
        if (!independent_set[u])
            continue;

        cost += g->W[u];
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (independent_set[v])
                return -1;
        }
    }
    return cost;
}

int main(int argc, char **argv)
{
    FILE *f = fopen(argv[1], "r");
    graph *g = graph_parse(f);
    fclose(f);

    printf("%d %d\n", g->N, g->V[g->N] / 2);

    long long offset = 0;
    int *A = malloc(sizeof(int) * g->N);
    int *S = malloc(sizeof(int) * g->N);
    int *reverse_mapping = malloc(sizeof(int) * g->N);

    for (int i = 0; i < g->N; i++)
    {
        A[i] = 1;
        S[i] = 0;
    }

    kernelize_csr(g->N, g->V, g->E, g->W, A, S, &offset, 7,
                  reduction_neighborhood_csr,
                  reduction_clique_csr,
                  reduction_domination_csr,
                  reduction_single_edge_csr,
                  reduction_extended_single_edge_csr,
                  //   reduction_twin_csr,
                  reduction_extended_twin_csr, // includes twin
                  reduction_unconfined_csr);

    int Nr = 0, Mr = 0;
    for (int u = 0; u < g->N; u++)
    {
        if (!A[u])
            continue;

        Nr++;
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (A[v])
                Mr++;
        }
    }

    printf("%d %d %lld\n", Nr, Mr / 2, offset);
    if (Nr != 0)
    {
        int run_pils = 16;
        graph *kernel = graph_subgraph(g, A, reverse_mapping);
        // FILE *f = fopen("kernel.graph", "w");
        // graph_store(f, kernel);
        // fclose(f);
        // graph *kernel = g; 
        // offset = 0;
        if (run_pils > 0)
        {
            pils *p = pils_init(kernel, run_pils);
            p->step_full = 1;
            p->step_reduced = 1;

            pils_run(kernel, p, 10000, 1, offset);
            long long w = mwis_validate(kernel, pils_get_best_independent_set(p));
            printf("%lld\n", w);

            pils_free(p);
        }
        else
        {
            local_search *ls = local_search_init(kernel, time(NULL));

            local_search_explore(kernel, ls, 10000, 1, offset);
            long long w = mwis_validate(kernel, ls->independent_set) + offset;
            printf("%lld\n", w);

            local_search_free(ls);
            graph_free(kernel);
        }
    }

    free(A);
    free(S);
    graph_free(g);

    return 0;
}