#include "reductions.h"

#include <stdlib.h>
#include <assert.h>

#define MAX_DEGREE 128

int reduction_heavy_vertex_csr(reduction_data *R, int N, const int *V, const int *E,
                               const long long *W, const int *A, int u, int *nRed, int *reducable)
{
    assert(A[u]);

    int degree = 0;
    for (int i = V[u]; i < V[u + 1] && degree <= MAX_DEGREE; i++)
        if (A[E[i]])
            degree++;

    if (degree > MAX_DEGREE)
        return 0;

    tiny_solver_solve_neighbourhood(R->solver, 1.0, W[u], u, V, E, W, A);

    if (!R->solver->time_limit_exceeded && W[u] >= R->solver->independent_set_weight)
    {
        *nRed = 1;
        reducable[0] = u;
        return 1;
    }

    return 0;
}
