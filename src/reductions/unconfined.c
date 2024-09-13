#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

int reduction_unconfined_csr(reduction_data *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u, int *nRed, int *reducable)
{
    assert(A[u]);

    reduction_data *rp = (reduction_data *)R;

    int n = 0, m = 0;

    assert(rp->Nb >= 4);

    int *S = rp->T[0], *NS = rp->T[1];
    int *S_B = rp->TB[0], *NS_B = rp->TB[1], *NSI_B = rp->TB[2];

    S[n++] = u;
    S_B[u] = 1;
    NSI_B[u] = 1;
    for (int i = V[u]; i < V[u + 1]; i++)
    {
        int v = E[i];
        if (!A[v])
            continue;
        NS[m++] = v;
        NS_B[v] = 1;
        NSI_B[v] = 1;
    }

    int res = 0, first = 1;
    while (m > 0)
    {
        int v = NS[--m];
        NS_B[v] = 0;
        if (S_B[v])
            continue;

        long long sw = 0;
        if (first)
            sw = W[u];
        else
            for (int i = V[v]; i < V[v + 1] && sw <= W[v]; i++)
                if (A[E[i]] && S_B[E[i]])
                    sw += W[E[i]];

        if (sw > W[v])
            continue;

        int xn = 0, x = -1;
        long long ds = 0;
        for (int i = V[v]; i < V[v + 1] && (sw + ds <= W[v] || xn <= 1); i++)
        {
            int w = E[i];
            if (A[w] && !NSI_B[w])
            {
                xn++;
                x = w;
                ds += W[w];
            }
        }

        if (sw + ds <= W[v]) // Can reduce u
        {
            res = 1;
            break;
        }
        else if (sw + ds > W[v] && xn == 1) // Extend S
        {
            first = 0;
            S[n++] = x;
            S_B[x] = 1;
            NSI_B[x] = 1;
            for (int i = V[x]; i < V[x + 1]; i++)
            {
                int w = E[i];
                if (!A[w] || NS_B[w])
                    continue;
                NS[m++] = w;
                NS_B[w] = 1;
                NSI_B[w] = 1;
            }
        }
    }

    while (m > 0)
    {
        int v = NS[--m];
        NS_B[v] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        int v = S[i];
        for (int j = V[v]; j < V[v + 1]; j++)
            NSI_B[E[j]] = 0;
        S_B[v] = 0;
        NSI_B[v] = 0;
    }

    if (res)
    {
        *nRed = 1;
        reducable[0] = u;
        return -1;
    }
    return 0;
}