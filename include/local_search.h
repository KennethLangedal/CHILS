#pragma once

#include "graph.h"

typedef struct
{
    int *C;
    int *IS;
} local_search;

local_search local_search_init(graph g);

void local_search_free(local_search ls);

void local_search_add_vertex(graph g, local_search ls, int u);

void local_search_remove_vertex(graph g, local_search ls, int u);

void local_search_greedy(graph g, local_search ls, int *order);

void local_search_k_one(graph g, local_search ls, int *order);

void local_search_k_c(graph g, local_search ls, int *order);