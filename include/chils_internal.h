#pragma once

#include "local_search.h"

typedef struct
{
    int p;
    double step_time;
    long long step_count;

    long long cost;
    int size;
    double time;

    local_search **LS, **LS_core;

    graph *d_core;
    int *FM, *RM, *A;
    long long *S1, *S2;
} chils;

chils *chils_init(graph *g, int p, unsigned int seed);

void chils_free(chils *c);

void chils_run(graph *g, chils *c, double tl, long long cl, int verbose);

void chils_set_solution(graph *g, chils *c, int i, const int *I);

int *chils_get_best_independent_set(chils *c);
