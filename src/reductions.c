#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

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
                   int *A, int *IS, long long *offset, int Nr, ...)
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
            IS[u] = 1;
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