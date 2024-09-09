#pragma once

#include <stdarg.h>

typedef int (*reduction_ptr)(void *R, int N, const int *V, const int *E,
                             const long long *W, const int *A, int u);

void kernelize_csr(int N, const int *V, const int *E, const long long *W,
                   int *A, int *IS, long long *offset, int Nr, ...);
