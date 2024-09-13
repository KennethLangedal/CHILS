#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

int reduction_neighborhood_csr(reduction_data *R, int N, const int *V, const int *E,
                               const long long *W, const int *A, int u, int *nRed, int *reducable)
{
    assert(A[u]);

    long long nw = 0;
    for (int i = V[u]; i < V[u + 1]; i++)
        if (A[E[i]])
            nw += W[E[i]];

    if (nw <= W[u])
    {
        *nRed = 1;
        reducable[0] = u;
        return 1;
    }
    return 0;
}
