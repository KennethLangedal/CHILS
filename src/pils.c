#include <omp.h>
#include <limits.h>
#include <stdlib.h>

#include "pils.h"
#include "reductions.h"

#define MIN_CORE 128

pils *pils_init(graph *g, int N)
{
    pils *p = malloc(sizeof(pils));

    p->N = N;

    p->step = 5.0;

    p->A = malloc(sizeof(int) * g->N);
    for (int i = 0; i < g->N; i++)
        p->A[i] = 0;

    p->LS = malloc(sizeof(local_search *) * N);

#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < N; i++)
            p->LS[i] = local_search_init(g, i);
    }

    return p;
}

void pils_free(pils *p)
{
    for (int i = 0; i < p->N; i++)
        local_search_free(p->LS[i]);
    free(p->LS);
    free(p->A);

    free(p);
}

void pils_print(pils *p, long long offset, double elapsed, int Nr)
{
    int best = 0, worst = 0;
    for (int i = 1; i < p->N; i++)
        if (p->LS[i]->cost > p->LS[best]->cost)
            best = i;
        else if (p->LS[i]->cost < p->LS[worst]->cost)
            worst = i;

    printf("\r%lld (%d) %lld (%d) %.2lf %d    ",
           p->LS[best]->cost + offset, best,
           p->LS[worst]->cost + offset, worst,
           elapsed, Nr);
    fflush(stdout);
}

void pils_run(graph *g, pils *p, double tl, int verbose, long long offset)
{
    double start = omp_get_wtime();
    double end = omp_get_wtime();
    double elapsed = end - start;

    if (verbose)
    {
        printf("Running PILS for %.2lf seconds\n", tl);
        pils_print(p, offset, elapsed, g->N);
    }

    graph *kernel = NULL;
    int *reverse_map = malloc(sizeof(int) * g->N);
    int *A = malloc(sizeof(int) * g->N);
    int Nr;

#pragma omp parallel shared(elapsed, kernel, reverse_map, A, Nr)
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
                pils_print(p, offset, elapsed, g->N);
            }
        }

        while (elapsed < tl)
        {
#pragma omp for
            for (int i = 0; i < p->N; i++)
                local_search_explore(g, p->LS[i], p->step, 0, 0);

#pragma omp for reduction(+ : Nr)
            for (int i = 0; i < g->N; i++)
            {
                int c = 0;
                for (int j = 0; j < p->N; j++)
                    c += p->LS[j]->independent_set[i];
                A[i] = c > 0 && c < p->N;
                Nr += A[i];
            }

#pragma omp single
            {
                end = omp_get_wtime();
                elapsed = end - start;
                if (verbose)
                    pils_print(p, offset, elapsed, g->N);
            }

            int best = 0, worst = 0;
            for (int i = 1; i < p->N; i++)
                if (p->LS[i]->cost > p->LS[best]->cost)
                    best = i;
                else if (p->LS[i]->cost < p->LS[worst]->cost)
                    worst = i;

            if (Nr > MIN_CORE)
            {

#pragma omp barrier
#pragma omp single
                {
                    kernel = graph_subgraph(g, A, reverse_map);
                }

#pragma omp for
                for (int i = 0; i < p->N; i++)
                {
                    local_search *ls_kernel = local_search_init(kernel, i);
                    long long ref = 0;
                    for (int u = 0; u < kernel->N; u++)
                        if (p->LS[i]->independent_set[reverse_map[u]])
                            ref += kernel->W[u];

                    local_search_explore(kernel, ls_kernel, p->step * 0.1, 0, 0);

                    if (i != best || ref <= ls_kernel->cost)
                        for (int u = 0; u < kernel->N; u++)
                            if (ls_kernel->independent_set[u] && !p->LS[i]->independent_set[reverse_map[u]])
                                local_search_add_vertex(g, p->LS[i], reverse_map[u]);

                    local_search_free(ls_kernel);
                }
#pragma omp single
                {
                    graph_free(kernel);
                    // end = omp_get_wtime();
                    // elapsed = end - start;
                }
            }

//             if (Nr <= MIN_CORE || p->LS[best]->cost == p->LS[worst]->cost)
//             {
// #pragma omp barrier
// #pragma omp for
//                 for (int i = 0; i < p->N; i++)
//                     if (i != best)
//                         local_search_scramble(g, p->LS[i], 128);
//             }

#pragma omp barrier
#pragma omp single
            {
                end = omp_get_wtime();
                elapsed = end - start;
                if (verbose)
                    pils_print(p, offset, elapsed, Nr);
                Nr = 0;
            }
        }

        //             int best = 0, worst = 0;
        //             for (int i = 1; i < p->N; i++)
        //                 if (p->LS[i]->cost > p->LS[best]->cost)
        //                     best = i;
        //                 else if (p->LS[i]->cost < p->LS[worst]->cost)
        //                     worst = i;

        //             if (kernel->N < MIN_CORE || p->LS[best]->cost == p->LS[worst]->cost)
        //             {

        // #pragma omp barrier
        // #pragma omp for
        //                 for (int i = 0; i < p->N; i++)
        //                 {
        //                     if (i == best)
        //                         continue;

        //                     p->LS[i]->log_enabled = 0;
        //                     for (int u = 0; u < g->N; u++)
        //                         if (p->LS[i]->independent_set[u])
        //                             local_search_remove_vertex(g, p->LS[i], u);
        //                     p->LS[i]->log_enabled = 1;
        //                 }
        //             }

        // #pragma omp barrier
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
