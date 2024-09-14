#pragma once

#include <stdarg.h>

// Buffers and bitvectors
typedef struct
{
    int Nb;
    int **T, **TB; // TB should be reset to 0 if used
} reduction_data;

typedef int (*reduction_ptr)(reduction_data *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u, int *Nr, int *reducable);

void kernelize_csr(int N, const int *V, const int *E, const long long *W,
                   int *A, int *S, long long *offset, int Nr, ...);

/**
 * Reduction rules for MWIS
 *
 * Each rule returns:
 *      1 if vertices can be included
 *      0 if the reduction failed
 *      -1 if vertices can be excluded
 **/

reduction_data *reduction_init(int N, int M);

void reduction_free(reduction_data *R);

int reduction_neighborhood_csr(reduction_data *R, int N, const int *V, const int *E,
                               const long long *W, const int *A, int u, int *nRed, int *reducable);

int reduction_clique_csr(reduction_data *R, int N, const int *V, const int *E,
                         const long long *W, const int *A, int u, int *nRed, int *reducable);

int reduction_domination_csr(reduction_data *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u, int *nRed, int *reducable);

int reduction_single_edge_csr(reduction_data *R, int N, const int *V, const int *E,
                         const long long *W, const int *A, int u, int* nRed, int* reducable);

int reduction_extended_single_edge_csr(reduction_data *R, int N, const int *V, const int *E,
                         const long long *W, const int *A, int u, int* nRed, int* reducable);

int reduction_unconfined_csr(reduction_data *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u, int *nRed, int *reducable);