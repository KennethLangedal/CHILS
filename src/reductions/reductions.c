#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <omp.h>

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
                   int *A, int *S, long long *offset, double tl, int Nr, ...)
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

    reduction_data *R = reduction_init(N, V[N]);
    int *reducable = malloc(sizeof(int) * N);
    int nRed = 0;
    double start = omp_get_wtime(), elapsed = 0.0;

    int rr = 0;
    while (rr < Nr && elapsed < tl)
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

        int res = reduction_rules[rr](R, N, V, E, W, A, u, &nRed, reducable);

        if (res == 1)
        {
            for (int i = 0; i < nRed; i++)
            {
                assert(A[reducable[i]]);
                int v = reducable[i];
                A[v] = 0;
                S[v] = 1;
                *offset += W[v];

                for (int i = V[v]; i < V[v + 1]; i++)
                {
                    A[E[i]] = 0;
                    kernelization_push_queue_dist_one(V, E, A, E[i], Nr, in_queue, queue, queue_size);
                }
            }
        }
        else if (res == -1)
        {
            for (int i = 0; i < nRed; i++)
            {
                assert(A[reducable[i]]);
                A[reducable[i]] = 0;
                kernelization_push_queue_dist_one(V, E, A, reducable[i], Nr, in_queue, queue, queue_size);
            }
        }
        nRed = 0;

        if (res != 0)
            rr = 0;

        elapsed = omp_get_wtime() - start;
    }

    free(reduction_rules);

    for (int i = 0; i < Nr; i++)
    {
        free(in_queue[i]);
        free(queue[i]);
    }

    free(in_queue);
    free(queue);
    free(queue_size);

    free(reducable);

    reduction_free(R);
}

reduction_data *reduction_init(int N, int M)
{
    reduction_data *rp = malloc(sizeof(reduction_data));

    rp->Nb = NUM_BUFFERS;
    rp->T = malloc(sizeof(int *) * rp->Nb);
    rp->TB = malloc(sizeof(int *) * rp->Nb);
    rp->solver = tiny_solver_init(N);

    for (int i = 0; i < rp->Nb; i++)
    {
        rp->T[i] = malloc(sizeof(int) * N);
        rp->TB[i] = malloc(sizeof(int) * N);

        for (int j = 0; j < N; j++)
            rp->TB[i][j] = 0;
    }

    return rp;
}

void reduction_free(reduction_data *R)
{
    reduction_data *rp = (reduction_data *)R;

    for (int i = 0; i < rp->Nb; i++)
    {
        free(rp->T[i]);
        free(rp->TB[i]);
    }

    free(rp->T);
    free(rp->TB);
    tiny_solver_free(rp->solver);

    free(rp);
}
