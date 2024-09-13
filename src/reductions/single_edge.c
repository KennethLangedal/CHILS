#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

int reduction_single_edge_csr(reduction_data *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u, int *nRed, int *reducable)
{
    assert(A[u]);

    reduction_data *rp = (reduction_data *)R;
    int *neighborhood_set = rp->TB[0];

    for (int i = V[u]; i < V[u + 1]; i++)
    {
        if (A[E[i]])
            neighborhood_set[E[i]] = 1;
    }

    for (int i = V[u]; i < V[u + 1]; i++)
    {
        int v = E[i];
        if (!A[v])
            continue;

        if (W[v] < W[u])
            continue;

        long long weight_sum = 0;
        for (int j = V[v]; j < V[v + 1]; j++)
        {
            if (A[E[j]] && !neighborhood_set[E[j]])
            {
                weight_sum += W[E[j]];
                if (weight_sum > W[v])
                    break;
            }
        }
        if (weight_sum <= W[v]) 
        {
            *nRed = 1;
            reducable[0] = u;

            // reset neighborhood set
            for (int i = V[u]; i < V[u + 1]; i++)
                neighborhood_set[E[i]] = 0;

            return -1;
        }
    }

    // reset neighborhood set
    for (int i = V[u]; i < V[u + 1]; i++)
        neighborhood_set[E[i]] = 0;
    return 0;
}
