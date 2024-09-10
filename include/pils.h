#pragma once

#include "local_search.h"

typedef struct
{
    int N;
    double step_full, step_reduced;

    int *C;
    local_search **LS;
} pils;

pils *pils_init(graph *g, int N);

void pils_free(pils *p);

void pils_run(graph *g, pils *p, double tl, int verbose);

int *pils_get_best_independent_set(pils *p);