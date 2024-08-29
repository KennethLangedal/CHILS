#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

typedef struct
{
    int *S, *NS, *T;
    int *S_B, *NSI_B, *IS_B;
} reduction_data;

void *reduction_init(int N, int M)
{
    reduction_data *rp = malloc(sizeof(reduction_data));

    rp->S = malloc(sizeof(int) * N);
    rp->NS = malloc(sizeof(int) * N);
    rp->T = malloc(sizeof(int) * N);
    rp->S_B = malloc(sizeof(int) * N);
    rp->NSI_B = malloc(sizeof(int) * N);
    rp->IS_B = malloc(sizeof(int) * N);

    for (int i = 0; i < N; i++)
    {
        rp->S_B[i] = 0;
        rp->NSI_B[i] = 0;
        rp->IS_B[i] = 0;
    }

    return rp;
}

void reduction_free(void *R)
{
    reduction_data *rp = (reduction_data *)R;

    free(rp->S);
    free(rp->NS);
    free(rp->T);
    free(rp->S_B);
    free(rp->NSI_B);
    free(rp->IS_B);
    free(rp);
}

int reduction_neighborhood_csr(void *R, int N, const int *V, const int *E,
                               const long long *W, const int *A, int u)
{
    if (!A[u])
        return 0;

    long long nw = 0;
    for (int i = V[u]; i < V[u + 1]; i++)
        if (A[E[i]])
            nw += W[E[i]];

    return nw <= W[u];
}

int reduction_unconfined_csr(void *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u)
{
    if (!A[u])
        return 0;

    reduction_data *rp = (reduction_data *)R;

    int n = 0, m = 0;

    rp->S[n++] = u;
    rp->S_B[u] = 1;
    rp->NSI_B[u] = 1;
    for (int i = V[u]; i < V[u + 1]; i++)
    {
        if (!A[E[i]])
            continue;
        rp->NS[m++] = E[i];
        rp->NSI_B[E[i]] = 1;
    }

    int res = 0, first = 1;
    while (m > 0)
    {
        int v = rp->NS[--m];
        if (rp->S_B[v])
            continue;

        long long sw = 0;
        if (first)
        {
            sw = W[u];
        }
        else
        {
            for (int i = V[v]; i < V[v + 1] && sw <= W[v]; i++)
                if (A[E[i]] && rp->S_B[E[i]])
                    sw += W[E[i]];
        }

        if (sw > W[v])
            continue;

        int p = 0;
        long long is = 0, min = LLONG_MAX;

        for (int i = V[v]; i < V[v + 1] && sw + is - min <= W[v]; i++)
        {
            int w = E[i];
            if (A[w] && !rp->NSI_B[w])
            {
                rp->T[p++] = w;
                is += W[w];
                rp->IS_B[w] = 1;
                if (W[w] < min)
                    min = W[w];
            }
        }

        int ind = sw + is - min <= W[v];
        for (int i = 0; i < p && ind; i++)
        {
            int w = rp->T[i];
            for (int j = V[w]; j < V[w + 1] && ind; j++)
                if (A[E[j]] && rp->IS_B[E[j]])
                    ind = 0;
        }

        for (int i = 0; i < p; i++)
            rp->IS_B[rp->T[i]] = 0;

        if (ind && sw + is <= W[v]) // Can reduce u
        {
            res = 1;
            break;
        }
        else if (ind && sw + is > W[v] && sw + is - min <= W[v]) // Extend S
        {
            first = 0;

            for (int i = 0; i < p; i++)
            {
                int w = rp->T[i];

                rp->S[n++] = w;
                rp->S_B[w] = 1;
                rp->NSI_B[w] = 1;
                for (int j = V[w]; j < V[w + 1]; j++)
                {
                    if (!A[E[j]])
                        continue;
                    rp->NS[m++] = E[j];
                    rp->NSI_B[E[j]] = 1;
                }
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        int v = rp->S[i];
        for (int j = V[v]; j < V[v + 1]; j++)
            rp->NSI_B[E[j]] = 0;
        rp->S_B[v] = 0;
        rp->NSI_B[v] = 0;
    }

    return res;
}