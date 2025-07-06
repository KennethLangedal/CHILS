#pragma once

#include "graph.h"

typedef struct
{
    // Solution
    long long cost, size;
    double time, time_ref;
    int *independent_set;

    // Queue structures
    int queue_count;
    int *queue, *in_queue;
    int *prev_queue, *in_prev_queue;

    // Graph structures
    int pool_size, max_queue;
    long long *adjacent_weight;
    int *tabu, *tightness, *temp, *mask, *pool;

    // Action log
    int log_count, log_alloc, log_enabled;
    int *log;

    unsigned int seed;
} local_search;

local_search *local_search_init(graph *g, unsigned int seed);

void local_search_free(local_search *ls);

void local_search_reset(graph *g, local_search *ls);

void local_search_in_order_solution(graph *g, local_search *ls);

void local_search_add_vertex(graph *g, local_search *ls, int u);

void local_search_remove_vertex(graph *g, local_search *ls, int u);

void local_search_lock_vertex(graph *g, local_search *ls, int u);

void local_search_unlock_vertex(graph *g, local_search *ls, int u);

void local_search_aap(graph *g, local_search *ls, int u, int imp);

void local_search_greedy(graph *g, local_search *ls);

void local_search_perturbe(graph *g, local_search *ls);

void local_search_explore(graph *g, local_search *ls, double tl, long long il, int verbose);

void local_search_unwind(graph *g, local_search *ls, int t);
