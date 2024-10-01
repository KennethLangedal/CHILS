#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

#define MAX_DEGREE 512
#define MAX_SOLVE 128

int reduction_unconfined_csr(reduction_data *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u, int *nRed, int *reducable)
{
    assert(A[u]);

    reduction_data *rp = (reduction_data *)R;

    int n = 0, m = 0;

    assert(rp->Nb >= 3);

    int *S = rp->T[0], *NS = rp->T[1];
    int *S_B = rp->TB[0], *NS_B = rp->TB[1], *NIS_B = rp->TB[2];

    S[n++] = u;
    S_B[u] = 1;
    NIS_B[u] = 1;
    for (int i = V[u]; i < V[u + 1] && m <= MAX_DEGREE; i++)
    {
        int v = E[i];
        if (!A[v])
            continue;
        NS[m++] = v;
        NS_B[v] = 1;
        NIS_B[v] = 1;
    }

    int res = 0, first = 1;
    while (m > 0 && m <= MAX_DEGREE)
    {
        int v = NS[--m];
        NS_B[v] = 0;
        if (S_B[v])
            continue;

        long long sw = 0, dw = 0;
        int dn = 0, x = -1;

        if (first && W[v] < W[u]) // Not a child
            continue;

        for (int i = V[v]; i < V[v + 1]; i++)
        {
            int w = E[i];
            if (!A[w])
                continue;

            if (S_B[w])
                sw += W[w];
            else if (!NIS_B[w])
            {
                dw += W[w];
                dn++;
                x = w;
            }
        }

        if (W[v] < sw) // Not a chld
            continue;

        if (sw + dw <= W[v]) // Condition 1. can reduce u
        {
            res = 1;
            break;
        }
        else if (sw + dw > W[v] && dn == 1) // Extending child
        {
            first = 0;
            S[n++] = x;
            S_B[x] = 1;
            NIS_B[x] = 1;
            for (int i = V[x]; i < V[x + 1]; i++)
            {
                int w = E[i];
                if (!A[w] || NS_B[w])
                    continue;
                NS[m++] = w;
                NS_B[w] = 1;
                NIS_B[w] = 1;
            }
        }
    }

    while (m > 0)
        NS_B[NS[--m]] = 0;

    for (int i = 0; i < n; i++)
    {
        int v = S[i];
        for (int j = V[v]; j < V[v + 1]; j++)
            NIS_B[E[j]] = 0;
        S_B[v] = 0;
        NIS_B[v] = 0;
    }

    if (res)
    {
        *nRed = 1;
        reducable[0] = u;
        return -1;
    }
    return 0;
}

int reduction_extended_unconfined_csr(reduction_data *R, int N, const int *V, const int *E,
                                      const long long *W, const int *A, int u, int *nRed, int *reducable)
{
    assert(A[u]);

    int n = 0, m = 0;

    assert(R->Nb >= 4);

    int *S = R->T[0], *NS = R->T[1], *T = R->T[2], *I = R->T[3];
    int *S_B = R->TB[0], *NS_B = R->TB[1], *NIS_B = R->TB[2];

    S[n++] = u;
    S_B[u] = 1;
    NIS_B[u] = 1;
    for (int i = V[u]; i < V[u + 1] && m <= MAX_DEGREE; i++)
    {
        int v = E[i];
        if (!A[v])
            continue;
        NS[m++] = v;
        NS_B[v] = 1;
        NIS_B[v] = 1;
    }

    int res = 0, first = 1;
    while (m > 0 && m <= MAX_DEGREE)
    {
        int v = NS[--m];
        NS_B[v] = 0;
        if (S_B[v])
            continue;

        long long sw = 0, dw = 0, dsw = 0;
        int dn = 0, in = 0;

        if (first && W[v] < W[u]) // Not a child
            continue;

        for (int i = V[v]; i < V[v + 1]; i++)
        {
            int w = E[i];
            if (!A[w])
                continue;

            if (S_B[w])
                sw += W[w];
            else if (!NIS_B[w])
            {
                dw += W[w];
                T[dn++] = w;
            }
        }

        if (W[v] < sw || dn > MAX_SOLVE) // Not a chld
            continue;

        if (W[v] >= sw + dw) // Can reduce u
        {
            res = 1;
            break;
        }

        tiny_solver_solve_subgraph(R->solver, 1.0, W[v] - sw, dn, T, V, E, W, A);
        if (R->solver->time_limit_exceeded)
            continue;

        dw = R->solver->independent_set_weight;
        dsw = 0;

        if (W[v] >= sw + dw) // Can reduce u
        {
            res = 1;
            break;
        }

        // Copy condidates
        for (int i = 0; i < dn; i++)
            I[i] = R->solver->independent_set[i] == 1 ? 1 : 0;

        // Find vertices to extend S, solve without each vertex included
        for (int i = 0; i < dn; i++)
        {
            if (!I[i])
                continue;

            int w = T[i];
            T[i] = T[dn - 1];
            tiny_solver_solve_subgraph(R->solver, 1.0, W[v] - sw, dn - 1, T, V, E, W, A);

            if (!R->solver->time_limit_exceeded && W[v] >= sw + R->solver->independent_set_weight)
            {
                first = 0;

                S[n++] = w;
                S_B[w] = 1;
                NIS_B[w] = 1;
                for (int j = V[w]; j < V[w + 1]; j++)
                {
                    int x = E[j];
                    if (!A[x] || NS_B[x])
                        continue;
                    NS[m++] = x;
                    NS_B[x] = 1;
                    NIS_B[x] = 1;
                }
            }

            T[i] = w;
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
            NIS_B[E[j]] = 0;
        S_B[v] = 0;
        NIS_B[v] = 0;
    }

    if (res)
    {
        *nRed = 1;
        reducable[0] = u;
        return -1;
    }

    return 0;
}