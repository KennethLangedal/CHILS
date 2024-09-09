#include "kernelization.h"
#include "reductions.h"

#include <stdlib.h>
#include <stdio.h>

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