#pragma once

#include "local_search.h"

typedef struct
{
    int N;
    double step;

    long long cost;
    double time;

    int *A;
    local_search **LS;
} chils;

chils *chils_init(graph *g, int N);

void chils_free(chils *p);

void chils_run(graph *g, chils *p, double tl, int verbose);

void chils_set_solution(graph *g, chils *p, const int *independent_set);

int *chils_get_best_independent_set(chils *p);
