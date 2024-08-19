#include "reductions.h"

#include <stdlib.h>

typedef struct
{
    int *S, *_S;
    int *NS_B, *S_B;
} reduction_data;

void *reduction_init(int N, int M)
{
    reduction_data *rp = malloc(sizeof(reduction_data));

    rp->S = malloc(sizeof(int) * N);
    rp->_S = malloc(sizeof(int) * N);
    rp->NS_B = malloc(sizeof(int) * N);
    rp->S_B = malloc(sizeof(int) * N);
    for (int i = 0; i < N; i++)
    {
        rp->NS_B[i] = 0;
        rp->S_B[i] = 0;
    }

    return rp;
}

void reduction_free(void *R)
{
    reduction_data *rp = (reduction_data *)R;

    free(rp->S);
    free(rp->_S);
    free(rp->NS_B);
    free(rp->S_B);
    free(rp);
}

int reduction_unconfined_csr(int N, const int *V, const int *E, const long long *W, const int *A, int u, void *R)
{
    reduction_data *rp = (reduction_data *)R;

    // int ns = 0;
    // for (int i = V[u]; i < V[u + 1]; i++)
    //     if (A[E[i]])
    //         rp->NS[ns++] = E[i];

    // for (int i = 0; i < ns; i++)
    // {
    //     int v = rp->NS[i];
    //     if (W[v] < W[u])
    //         continue;
    // }

    return 1;
}