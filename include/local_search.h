#pragma once

#include "graph.h"

#define ADD 1
#define REMOVE 0

#define MAX_LOG 100000000

typedef struct
{
    // Cost and solution
    int c;
    int *IS;

    // Search structures
    int qc, a;
    int *NW, *T, *Q, *C, *_Q, *_C; // TODO, add TABU

    // Action log
    int lc;
    int *A, *L;
} local_search;

local_search *local_search_init(graph g);

void local_search_free(local_search *ls);

void local_search_add_vertex(graph g, local_search *ls, int u);

void local_search_remove_vertex(graph g, local_search *ls, int u);

void local_search_greedy(graph g, local_search *ls);

void local_search_unwind(graph g, local_search *ls, int t);

void local_search_shuffle_queue(local_search *ls);
