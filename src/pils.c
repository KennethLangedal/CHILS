#include <omp.h>
#include <stdlib.h>

#include "pils.h"
#include "reductions.h"

pils *pils_init(graph *g, int N)
{
    pils *p = malloc(sizeof(pils));

    p->N = N;

    p->step = 5.0;

    p->C = malloc(sizeof(int) * g->N);
    p->A = malloc(sizeof(int) * g->N);
    p->S = malloc(sizeof(int) * g->N);

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
    free(p->C);
    free(p->A);
    free(p->S);

    free(p);
}

void pils_run(graph *g, pils *p, double tl, int verbose, long long offset)
{
    double start, end;
    start = omp_get_wtime();

    if (verbose)
    {
        printf("Running PILS for %.2lf seconds\n", tl);

        long long best = 0;
        for (int i = 0; i < p->N; i++)
            if (p->LS[i]->cost > best)
                best = p->LS[i]->cost;

        printf("\r%lld %.2lf 0    ", best + offset, 0.0);
        fflush(stdout);
    }

    end = omp_get_wtime();
    double elapsed = end - start;

#pragma omp parallel shared(elapsed)
    {
#pragma omp for
        for (int i = 0; i < p->N; i++)
        {
            if (p->LS[i]->cost == 0 && i == 0)
                local_search_in_order_solution(g, p->LS[i]);

            local_search_greedy(g, p->LS[i]);
        }

#pragma omp single
        {
            end = omp_get_wtime();
            elapsed = end - start;
        }

        if (verbose)
        {
#pragma omp single
            {
                long long best = 0;
                for (int i = 0; i < p->N; i++)
                    if (p->LS[i]->cost > best)
                        best = p->LS[i]->cost;

                printf("\r%lld %.2lf %d    ", best + offset, elapsed, g->N);
                fflush(stdout);
            }
        }

        while (elapsed < tl)
        {
#pragma omp for
            for (int i = 0; i < g->N; i++)
                p->C[i] = 0;

#pragma omp for
            for (int i = 0; i < g->N; i++)
                for (int j = 0; j < p->N; j++)
                    p->C[i] += p->LS[j]->independent_set[i];

            int Nr = 0;
#pragma omp for
            for (int i = 0; i < p->N; i++)
            {
                Nr = 0;
                for (int j = 0; j < g->N; j++)
                {
                    if (p->C[j] == 0 || p->C[j] == p->N)
                        p->LS[i]->tabu[j] = 1;
                    else
                        Nr++;
                }
            }

#pragma omp for
            for (int i = 0; i < p->N; i++)
                local_search_explore(g, p->LS[i], p->step * ((double)Nr / (double)g->N), 0, 0);

#pragma omp for
            for (int i = 0; i < p->N; i++)
                for (int j = 0; j < g->N; j++)
                    p->LS[i]->tabu[j] = 0;

#pragma omp for
            for (int i = 0; i < p->N; i++)
                local_search_explore(g, p->LS[i], p->step, 0, 0);

#pragma omp single
            {
                end = omp_get_wtime();
                elapsed = end - start;
            }

            if (verbose)
            {
#pragma omp single
                {
                    long long best = 0;
                    for (int i = 0; i < p->N; i++)
                        if (p->LS[i]->cost > best)
                            best = p->LS[i]->cost;

                    printf("\r%lld %.2lf %d    ", best + offset, elapsed, Nr);
                    fflush(stdout);
                }
            }
        }
    }

    if (verbose)
        printf("\n");
}

void pils_set_solution(graph *g, pils *p, const int *independent_set)
{
#pragma omp for
    for (int i = 0; i < p->N; i++)
        for (int j = 0; j < g->N; j++)
            if (independent_set[j])
                local_search_add_vertex(g, p->LS[i], j);
}

int *pils_get_best_independent_set(pils *p)
{
    int best = 0;
    for (int i = 0; i < p->N; i++)
        if (p->LS[i]->cost > p->LS[best]->cost)
            best = i;

    return p->LS[best]->independent_set;
}
