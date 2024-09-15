#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

int reduction_twin_csr(reduction_data *R, int N, const int *V, const int *E,
                       const long long *W, const int *A, int u, int *nRed, int *reducable)
{
    assert(A[u]);

    reduction_data *rp = (reduction_data *)R;
    int *neighborhood_set = rp->TB[0];
    int *candidates = rp->T[0];
    long long neighborhood_weight = 0;
    int min_deg_neighbor = -1;
    int min_deg = N;
    int deg_u = 0;
    int n_candidates = 0;

    for (int i = V[u]; i < V[u + 1]; i++)
    {
        if (A[E[i]])
        {
            neighborhood_set[E[i]] = 1;
            neighborhood_weight += W[E[i]];
            deg_u++;
            if (min_deg > V[E[i] + 1] - V[E[i]])
            {
                min_deg = V[E[i] + 1] - V[E[i]];
                min_deg_neighbor = E[i];
            }
        }
    }

    // check the neighborhood of min degree neighbor for twin candidates 
    for (int i = V[min_deg_neighbor]; i < V[min_deg_neighbor + 1]; i++)
    {
        if (A[E[i]] && E[i] != u && !neighborhood_set[E[i]] && V[E[i] + 1] - V[E[i]] >= deg_u)
            candidates[n_candidates++] = E[i];
    }

    // check if twin candidates are actual twins i.e. have the same neighborhood
    for (int i = 0; i < n_candidates; i++)
    {
        int twin_neighbor_count = 0;
        int twin = candidates[i];
        assert(A[twin]);
        int is_twin = 1;
        for (int j = V[twin]; j < V[twin + 1]; j++)
        {
            if (!A[E[j]])
                continue;
            twin_neighbor_count++;

            if(!neighborhood_set[E[j]] || twin_neighbor_count > deg_u)
            {
                is_twin = 0;
                break; 
            }
        }
        if (!is_twin || twin_neighbor_count != deg_u)
            continue;
        
        if (W[twin] + W[u] >= neighborhood_weight) 
        {
            *nRed = 2;
            reducable[0] = u;
            reducable[1] = twin;

            for (int i = V[u]; i < V[u + 1]; i++)
                neighborhood_set[E[i]] = 0;

            return 1;
        }
        // TODO: else fold
    }

    // reset neighborhood set
    for (int i = V[u]; i < V[u + 1]; i++)
        neighborhood_set[E[i]] = 0;
    return 0;
}
