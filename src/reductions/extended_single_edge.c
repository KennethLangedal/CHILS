#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

int reduction_extended_single_edge_csr(reduction_data *R, int N, const int *V, const int *E,
                                       const long long *W, const int *A, int u, int *nRed, int *reducable)
{
    assert(A[u]);

    reduction_data *rp = (reduction_data *)R;
    int *neighborhood_set = rp->TB[0];
    long long neighborhood_weight = 0;

    for (int i = V[u]; i < V[u + 1]; i++)
    {
        if (A[E[i]])
        {
            neighborhood_set[E[i]] = 1;
            neighborhood_weight += W[E[i]];
        }
    }

    int n = 0;
    for (int i = V[u]; i < V[u + 1]; i++)
    {
        int v = E[i];
        if (!A[v])
            continue;

        if (W[u] < neighborhood_weight - W[v])
            continue;

        for (int j = V[v]; j < V[v + 1]; j++)
        {
            if (!A[E[j]] || u == E[j])
                continue;
            if (neighborhood_set[E[j]])
            {
                reducable[n++] = E[j];
                neighborhood_set[E[j]] = 0; 
                neighborhood_weight -= W[E[j]];
            }
        }
    }

    // reset neighborhood set
    for (int i = V[u]; i < V[u + 1]; i++)
        neighborhood_set[E[i]] = 0;

    if (n > 0)
    {
        *nRed = n;
        return -1;
    }
    return 0;
}
