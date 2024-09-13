#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

int reduction_clique_csr(reduction_data *R, int N, const int *V, const int *E,
                         const long long *W, const int *A, int u, int *nRed, int *reducable)
{
    assert(A[u]);
    reduction_data *rp = (reduction_data *)R;
    int *isolated = rp->T[0];
    int *to_reduce = rp->T[1];
    int *clique_set = rp->TB[0];
    isolated[0] = u;
    int isolated_count = 1;

    int u_clique_count = 1;
    clique_set[u] = 1;
    for (int i = V[u]; i < V[u + 1]; i++)
    {
        if (A[E[i]])
        {
            clique_set[E[i]] = 1;
            u_clique_count++;
        }
    }

    int is_clique = 1;
    for (int i = V[u]; i < V[u + 1]; i++)
    {
        int c = E[i];
        if (!A[c])
            continue;

        int is_isolated = 1;
        int count = 1;
        for (int j = V[c]; j < V[c + 1]; j++)
        {
            if (!A[E[j]])
                continue;
            if (clique_set[E[j]])
                count++;
            else
                is_isolated = 0;
        }

        if (count != u_clique_count)
        {
            is_clique = 0;
            break;
        }

        // c is also isolated if it has the same number of neighbors as u
        if (is_isolated)
            isolated[isolated_count++] = c;
    }

    // reset clique set
    clique_set[u] = 0;
    for (int i = V[u]; i < V[u + 1]; i++)
        clique_set[E[i]] = 0;

    if (!is_clique)
        return 0;

    long long max_weight_isolated = 0;
    int max_weight_isolated_node = -1;
    for (int i = 0; i < isolated_count; i++)
    {
        if (W[isolated[i]] > max_weight_isolated)
        {
            max_weight_isolated = W[isolated[i]];
            max_weight_isolated_node = isolated[i];
        }
    }

    // get max weight of clique
    long long max_weight_clique = 0;
    int n_reduceable = 0;
    for (int i = V[u]; i < V[u + 1]; i++)
    {
        if (!A[E[i]])
            continue;

        if (W[E[i]] > max_weight_clique)
            max_weight_clique = W[E[i]];

        // reduce clique vertices of less weight than max_weight_isolated
        if (W[E[i]] <= max_weight_isolated && E[i] != max_weight_isolated_node)
            to_reduce[n_reduceable++] = E[i];
    }

    if (max_weight_isolated == max_weight_clique)
    {
        // include max_weight_isolated_node
        *nRed = 1;
        reducable[0] = max_weight_isolated_node;
        return 1;
    }
    else if (n_reduceable == 0)
    {
        return 0;
    }
    else
    {
        // exclude all vertices of the clique with weight less than max_weight_isolated
        *nRed = n_reduceable;
        for (int i = 0; i < n_reduceable; i++)
            reducable[i] = to_reduce[i];

        return -1;
    }
}
