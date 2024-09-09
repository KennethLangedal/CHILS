#pragma once

#include "clique_graph.h"

typedef struct
{
    // Solution
    long long cost;
    int *independent_set;

    // Queue structures
    int queue_count;
    int *queue, *in_queue;
    int *prev_queue, *in_prev_queue;

    // Graph structures
    int fsc;
    int *clique_status, *fast_set;

    // Action log
    int log_count;
    int *log, *action;

    unsigned int seed;
} clique_local_search;

clique_local_search *clique_local_search_init(clique_graph *g, unsigned int seed);

void clique_local_search_free(clique_local_search *ls);

void clique_local_search_in_order_solution(clique_graph *g, clique_local_search *ls);

void clique_local_search_add_vertex(clique_graph *g, clique_local_search *ls, int u);

void clique_local_search_remove_vertex(clique_graph *g, clique_local_search *ls, int u);

void clique_local_search_greedy(clique_graph *g, clique_local_search *ls);

void clique_local_search_explore(clique_graph *g, clique_local_search *ls, int it, int k, int verbose);

void clique_local_search_unwind(clique_graph *g, clique_local_search *ls, int t);
