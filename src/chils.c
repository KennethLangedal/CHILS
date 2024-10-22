#include <omp.h>
#include <limits.h>
#include <stdlib.h>

#include "chils.h"

#define MIN_CORE 512

chils *chils_init(graph *g, int N)
{
    chils *p = malloc(sizeof(chils));

    p->N = N;
    p->step = 5.0;

    p->cost = 0;
    p->time = 0.0;

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

void chils_free(chils *p)
{
    for (int i = 0; i < p->N; i++)
        local_search_free(p->LS[i]);
    free(p->LS);
    free(p->A);

    free(p);
}

void chils_print(chils *p, double elapsed, int Nr)
{
    int best = 0, worst = 0;
    for (int i = 1; i < p->N; i++)
        if (p->LS[i]->cost >= p->LS[best]->cost)
            best = i;
        else if (p->LS[i]->cost < p->LS[worst]->cost)
            worst = i;

    printf("\r%lld (%d %.2lf) %lld (%d %.2lf) %.2lf %d    ",
           p->LS[best]->cost, best, p->LS[best]->time,
           p->LS[worst]->cost, worst, p->LS[worst]->time,
           elapsed, Nr);
    fflush(stdout);
}

int chils_find_last_best(chils *p)
{
    int best = 0;
    for (int i = 1; i < p->N; i++)
        if (p->LS[i]->cost >= p->LS[best]->cost)
            best = i;
    return best;
}

void chils_update_best(chils *p)
{
    for (int i = 0; i < p->N; i++)
    {
        if (p->LS[i]->cost > p->cost ||
            (p->LS[i]->cost == p->cost && p->LS[i]->time < p->time))
        {
            p->cost = p->LS[i]->cost;
            p->time = p->LS[i]->time;
        }
    }
}

void chils_run(graph *g, chils *p, double tl, int verbose)
{
    double start = omp_get_wtime();
    double end = omp_get_wtime();
    double elapsed = end - start;

    if (verbose)
    {
        printf("Running chils for %.2lf seconds\n", tl);
        chils_print(p, elapsed, g->N);
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
            else if (p->LS[i]->cost == 0)
                local_search_add_vertex(g, p->LS[i], rand_r(&p->LS[i]->seed) % g->N);
        }

#pragma omp single
        {
            end = omp_get_wtime();
            elapsed = end - start;

            chils_update_best(p);

            if (verbose)
                chils_print(p, elapsed, g->N);
        }

        while (elapsed < tl)
        {
            /* Full graph LS */
#pragma omp for
            for (int i = 0; i < p->N; i++)
            {
                double remaining_time = tl - (omp_get_wtime() - start);
                double duration = p->step;
                if (remaining_time < duration)
                    duration = remaining_time;
                if (duration > 0.0)
                    local_search_explore(g, p->LS[i], duration, 0);
            }

            /* Mark the CHILS core */
#pragma omp for reduction(+ : Nr)
            for (int i = 0; i < g->N; i++)
            {
                int c = 0;
                for (int j = 0; j < p->N; j++)
                    c += p->LS[j]->independent_set[i];
                A[i] = c > 0 && c < p->N;
                Nr += A[i];
            }

            /* Find the best solution */
            int best = chils_find_last_best(p);

            /* Construct the CHILS core */
#pragma omp single
            {
                chils_update_best(p);
                kernel = graph_subgraph(g, A, reverse_map);
            }

            /* CHILS core LS */
#pragma omp for
            for (int i = 0; i < p->N; i++)
            {
                if (kernel->N == 0)
                    continue;

                double remaining_time = tl - (omp_get_wtime() - start);
                double duration = p->step * 0.5;
                if (remaining_time < duration)
                    duration = remaining_time;

                if (duration < 0.0)
                    continue;

                local_search *ls_kernel = local_search_init(kernel, i);

                long long ref = 0;
                for (int u = 0; u < kernel->N; u++)
                    if (p->LS[i]->independent_set[reverse_map[u]])
                        ref += kernel->W[u];

                local_search_explore(kernel, ls_kernel, duration, 0);

                if (ref <= ls_kernel->cost || (i != best && (i % 2) == 0))
                    for (int u = 0; u < kernel->N; u++)
                        if (ls_kernel->independent_set[u] && !p->LS[i]->independent_set[reverse_map[u]])
                            local_search_add_vertex(g, p->LS[i], reverse_map[u]);

                if (ref < ls_kernel->cost)
                    p->LS[i]->time = (omp_get_wtime() - p->LS[i]->time_ref) - (duration - ls_kernel->time);

                local_search_free(ls_kernel);
            }

            /* Find the best solution after LS on the CHILS core */
            best = chils_find_last_best(p);

#pragma omp single
            {
                chils_update_best(p);
                graph_free(kernel);
            }

#pragma omp for
            for (int i = 0; i < p->N; i++)
            {
                if (Nr < MIN_CORE && i != best && (i % 2) == 0)
                    local_search_perturbe(g, p->LS[i]);
            }

#pragma omp single
            {
                end = omp_get_wtime();
                elapsed = end - start;
                if (verbose)
                    chils_print(p, elapsed, Nr);
                Nr = 0;
            }
        }
    }

    if (verbose)
        printf("\n");

    free(reverse_map);
    free(A);
}

void chils_set_solution(graph *g, chils *p, const int *independent_set)
{
#pragma omp for
    for (int i = 0; i < p->N; i++)
        for (int j = 0; j < g->N; j++)
            if (independent_set[j])
                local_search_add_vertex(g, p->LS[i], j);
}

int *chils_get_best_independent_set(chils *p)
{
    int best = 0;
    for (int i = 0; i < p->N; i++)
        if (p->LS[i]->cost > p->LS[best]->cost)
            best = i;

    return p->LS[best]->independent_set;
}
