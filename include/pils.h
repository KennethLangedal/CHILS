#pragma once

#include "local_search.h"

typedef struct
{
    int N;
    double step;

    int *C, *A, *S;
    local_search **LS;
} pils;

pils *pils_init(graph *g, int N);

void pils_free(pils *p);

void pils_run(graph *g, pils *p, double tl, int verbose, long long offset);

void pils_set_solution(graph *g, pils *p, const int *independent_set);

int *pils_get_best_independent_set(pils *p);