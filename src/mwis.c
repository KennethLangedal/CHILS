#include "mwis.h"

long long mwis_validate(graph *g, int *independent_set)
{
    long long C = 0;
    for (int u = 0; u < g->N; u++)
    {
        if (!independent_set[u])
            continue;

        C += g->W[u];
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (independent_set[v])
                return -1;
        }
    }
    return C;
}