#include <omp.h>
#include <stdlib.h>

#include "pils.h"

pils *pils_init(graph *g, int N)
{
    pils *p = malloc(sizeof(pils));

    p->N = N;

    p->step_full = 5.0;
    p->step_reduced = 5.0;

    p->C = malloc(sizeof(int) * g->N);
    for (int i = 0; i < g->N; i++)
        p->C[i] = 0;

    p->LS = malloc(sizeof(local_search *) * N);

#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < N; i++)
        {
            p->LS[i] = local_search_init(g, i);
        }
    }

    return p;
}

void pils_free(pils *p)
{
    for (int i = 0; i < p->N; i++)
        local_search_free(p->LS[i]);
    free(p->LS);

    free(p);
}

void pils_run(graph *g, pils *p, double tl, int verbose, int offset)
{
    double start, end;
    start = omp_get_wtime();

    if (verbose)
    {
        printf("Running PILS using for %.2lf seconds\n", tl);

        long long best = 0;
        for (int i = 0; i < p->N; i++)
            if (p->LS[i]->cost > best)
                best = p->LS[i]->cost;
        printf("\r%lld %.2lf   ", best + offset, 0.0);
        fflush(stdout);
    }

#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < p->N; i++)
        {
            if (p->LS[i]->cost == 0 && i == 0)
                local_search_in_order_solution(g, p->LS[i]);

            local_search_greedy(g, p->LS[i]);
        }

        end = omp_get_wtime();
        double elapsed = end - start;

        while (elapsed < tl)
        {
#pragma omp for
            for (int i = 0; i < g->N; i++)
                p->C[i] = 0;

#pragma omp for
            for (int i = 0; i < g->N; i++)
                for (int j = 0; j < p->N; j++)
                    p->C[i] += p->LS[j]->independent_set[i];

#pragma omp for
            for (int i = 0; i < g->N; i++)
            {
                if (p->C[i] < p->N)
                    continue;

                for (int j = 0; j < p->N; j++)
                    local_search_lock_vertex(g, p->LS[j], i);
            }

#pragma omp for
            for (int i = 0; i < p->N; i++)
                local_search_explore(g, p->LS[i], p->step_reduced, 0, offset);

#pragma omp for
            for (int i = 0; i < g->N; i++)
            {
                if (p->C[i] < p->N)
                    continue;

                for (int j = 0; j < p->N; j++)
                    local_search_unlock_vertex(g, p->LS[j], i);
            }

#pragma omp for
            for (int i = 0; i < p->N; i++)
                local_search_explore(g, p->LS[i], p->step_full, 0, offset);

            end = omp_get_wtime();
            elapsed = end - start;
            if (verbose)
            {
#pragma omp single
                {
                    long long best = 0;
                    for (int i = 0; i < p->N; i++)
                        if (p->LS[i]->cost > best)
                            best = p->LS[i]->cost;

                    printf("\r%lld %.2lf   ", best + offset, elapsed);
                    fflush(stdout);
                }
            }
        }
    }

    if (verbose)
        printf("\n");
}

int *pils_get_best_independent_set(pils *p)
{
    int best = 0;
    for (int i = 0; i < p->N; i++)
        if (p->LS[i]->cost > p->LS[best]->cost)
            best = i;

    return p->LS[best]->independent_set;
}