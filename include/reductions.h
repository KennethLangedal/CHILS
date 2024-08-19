#pragma once

void *reduction_init(int N, int M);

void reduction_free(void *R);

int reduction_unconfined_csr(int N, const int *V, const int *E, const long long *W, const int *A, int u, void *R);