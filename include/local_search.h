#pragma once

#include "graph.h"

#define REMOVE 0
#define ADD 1
#define REMOVE_TABU 2
#define ADD_TABU 3

#define MAX_LOG 1000000

typedef struct
{
    // Cost and solution
    long long c;
    int *IS;

    // Search structures
    int qc;
    int *NW, *T, *P, *Q, *C, *_Q, *_C;

    // Action log
    int lc;
    int *A, *L;
} local_search;

local_search *local_search_init(graph g);

void local_search_free(local_search *ls);

void local_search_add_vertex(graph g, local_search *ls, int u);

void local_search_remove_vertex(graph g, local_search *ls, int u);

void local_search_add_vertex_tabu(graph g, local_search *ls, int u);

void local_search_remove_vertex_tabu(graph g, local_search *ls, int u);

void local_search_greedy(graph g, local_search *ls);

void local_search_unwind(graph g, local_search *ls, int t);

void local_search_unwind_tabu(graph g, local_search *ls, int t);
