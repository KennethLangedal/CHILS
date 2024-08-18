#pragma once

#include "graph.h"

#define REMOVE 0
#define ADD 1

#define MAX_LOG 10000000

typedef struct
{
    // Cost and solution
    long long c;
    int *IS;

    // Search structures
    int qc;
    int *NW, *T, *_T, *P, *_P, *Q, *C, *_Q, *_C;

    // Action log
    int lc;
    int *A, *L;

    unsigned int seed;
} local_search;

local_search *local_search_init(graph g, unsigned int seed);

void local_search_free(local_search *ls);

void local_search_add_vertex(graph g, local_search *ls, int u);

void local_search_remove_vertex(graph g, local_search *ls, int u);

void local_search_lock_vertex(graph g, local_search *ls, int u);

void local_search_unlock_vertex(graph g, local_search *ls, int u);

void local_search_aap(graph g, local_search *ls, int u);

void local_search_greedy(graph g, local_search *ls);

void local_search_unwind(graph g, local_search *ls, int t);
