#pragma once

#include "graph.h"

typedef struct
{
    int *C;
    int *IS, *nw, *tabu;
} local_search;

local_search local_search_init(graph g);

void local_search_copy(graph g, local_search src, local_search dest);

void local_search_free(local_search ls);

void local_search_add_vertex(graph g, local_search ls, int u);

void local_search_remove_vertex(graph g, local_search ls, int u);

void local_search_lock_in_vertex(graph g, local_search ls, int u);

void local_search_lock_out_vertex(graph g, local_search ls, int u);

void local_search_greedy(graph g, local_search ls, int *order);

void local_search_k_one(graph g, local_search ls, int *order);

void local_search_k_c(graph g, local_search ls, int *order);