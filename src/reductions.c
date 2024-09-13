#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

#define NUM_BUFFERS 4

static inline void kernelization_push_queue_dist_one(const int *V, const int *E, const int *A, int u,
                                                     int Nr, int **in_queue, int **queue, int *queue_size)
{
    for (int i = V[u]; i < V[u + 1]; i++)
    {
        int v = E[i];
        if (!A[v])
            continue;

        for (int k = 0; k < Nr; k++)
        {
            if (!in_queue[k][v])
            {
                in_queue[k][v] = 1;
                queue[k][queue_size[k]++] = v;
            }
        }
    }
}

void kernelize_csr(int N, const int *V, const int *E, const long long *W,
                   int *A, int *S, long long *offset, int Nr, ...)
{
    va_list argptr;
    va_start(argptr, Nr);

    reduction_ptr *reduction_rules = malloc(sizeof(reduction_ptr) * Nr);
    for (int i = 0; i < Nr; i++)
        reduction_rules[i] = va_arg(argptr, reduction_ptr);

    va_end(argptr);

    int **in_queue = malloc(sizeof(int *) * Nr);
    int **queue = malloc(sizeof(int *) * Nr);
    int *queue_size = malloc(sizeof(int) * Nr);

    for (int i = 0; i < Nr; i++)
    {
        in_queue[i] = malloc(sizeof(int) * N);
        queue[i] = malloc(sizeof(int) * N);
        queue_size[i] = N;

        for (int j = 0; j < N; j++)
        {
            in_queue[i][j] = 1;
            queue[i][j] = j;
        }
    }

    void *R = reduction_init(N, V[N]);

    int rr = 0;
    while (rr < Nr)
    {
        if (queue_size[rr] == 0)
        {
            rr++;
            continue;
        }

        int u = queue[rr][--queue_size[rr]];
        in_queue[rr][u] = 0;

        if (!A[u])
            continue;

        int res = reduction_rules[rr](R, N, V, E, W, A, u);

        if (res == 1)
        {
            A[u] = 0;
            S[u] = 1;
            *offset += W[u];

            for (int i = V[u]; i < V[u + 1]; i++)
                A[E[i]] = 0;

            for (int i = V[u]; i < V[u + 1]; i++)
                kernelization_push_queue_dist_one(V, E, A, E[i], Nr, in_queue, queue, queue_size);

            printf("\r%d  ", queue_size[1]);
            fflush(stdout);
        }
        else if (res == -1)
        {
            A[u] = 0;
            kernelization_push_queue_dist_one(V, E, A, u, Nr, in_queue, queue, queue_size);

            printf("\r%d  ", queue_size[1]);
            fflush(stdout);
        }

        if (res != 0)
            rr = 0;
    }
    printf("\n");

    free(reduction_rules);

    for (int i = 0; i < Nr; i++)
    {
        free(in_queue[i]);
        free(queue[i]);
    }

    free(in_queue);
    free(queue);
    free(queue_size);

    reduction_free(R);
}

// Buffers and bitvectors (should always be reset by reduction rule)
typedef struct
{
    int Nb;
    int **T, **TB;
} reduction_data;

void *reduction_init(int N, int M)
{
    reduction_data *rp = malloc(sizeof(reduction_data));

    rp->Nb = NUM_BUFFERS;
    rp->T = malloc(sizeof(int *) * rp->Nb);
    rp->TB = malloc(sizeof(int *) * rp->Nb);

    for (int i = 0; i < rp->Nb; i++)
    {
        rp->T[i] = malloc(sizeof(int) * N);
        rp->TB[i] = malloc(sizeof(int) * N);

        for (int j = 0; j < N; j++)
            rp->TB[i][j] = 0;
    }

    return rp;
}

void reduction_free(void *R)
{
    reduction_data *rp = (reduction_data *)R;

    for (int i = 0; i < rp->Nb; i++)
    {
        free(rp->T[i]);
        free(rp->TB[i]);
    }

    free(rp->T);
    free(rp->TB);

    free(rp);
}

int reduction_neighborhood_csr(void *R, int N, const int *V, const int *E,
                               const long long *W, const int *A, int u)
{
    assert(A[u]);

    long long nw = 0;
    for (int i = V[u]; i < V[u + 1]; i++)
        if (A[E[i]])
            nw += W[E[i]];

    if (nw <= W[u])
        return 1;
    return 0;
}

int reduction_unconfined_csr(void *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u)
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
        return -1;
    return 0;
}