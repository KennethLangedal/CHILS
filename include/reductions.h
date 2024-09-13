#pragma once

#include <stdarg.h>

typedef int (*reduction_ptr)(void *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u);

void kernelize_csr(int N, const int *V, const int *E, const long long *W,
                   int *A, int *S, long long *offset, int Nr, ...);

/**
 * Reduction rules for MWIS
 *
 * Each rule returns:
 *      1 if u can be included
 *      0 if the reduction failed
 *      -1 if u can be excluded
 **/

void *reduction_init(int N, int M);

void reduction_free(void *R);

int reduction_neighborhood_csr(void *R, int N, const int *V, const int *E,
                               const long long *W, const int *A, int u);

int reduction_unconfined_csr(void *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u);