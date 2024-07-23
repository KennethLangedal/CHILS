#pragma once

#include "graph.h"

typedef struct
{
    unsigned int N, *C;
    unsigned int *IS, *W, *V, *E;
} local_search_avx;

local_search_avx local_search_avx_init(graph g);

void local_search_avx_free(local_search_avx ls);

void local_search_avx_add_vertex(local_search_avx ls, int u);

void local_search_avx_remove_vertex(local_search_avx ls, int u);

void local_search_avx_greedy(local_search_avx ls, int *order);