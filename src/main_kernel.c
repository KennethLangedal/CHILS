#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "reductions.h"

int main(int argc, char **argv)
{
    FILE *f = fopen(argv[1], "r");
    graph *g = graph_parse(f);
    fclose(f);

    printf("%d %d\n", g->N, g->V[g->N] / 2);

    long long offset = 0;
    int *A = malloc(sizeof(int) * g->N);
    int *S = malloc(sizeof(int) * g->N);

    for (int i = 0; i < g->N; i++)
    {
        A[i] = 1;
        S[i] = 0;
    }

    kernelize_csr(g->N, g->V, g->E, g->W, A, S, &offset, 4,
                  reduction_neighborhood_csr,
                  reduction_clique_csr,
                  reduction_domination_csr,
                  reduction_unconfined_csr
                  );

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

    graph_free(g);

    return 0;
}