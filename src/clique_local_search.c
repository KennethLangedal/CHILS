#include "clique_local_search.h"

#include <stdlib.h>
#include <limits.h>
#include <assert.h>

#define REMOVE 0
#define ADD 1

#define MAX_LOG 10000000

clique_local_search *clique_local_search_init(clique_graph *g, unsigned int seed)
{
    clique_local_search *ls = malloc(sizeof(clique_local_search));

    ls->cost = 0;
    ls->independent_set = malloc(sizeof(int) * g->N);

    ls->queue_count = g->NC;
    ls->queue = malloc(sizeof(int) * g->NC);
    ls->in_queue = malloc(sizeof(int) * g->NC);
    ls->prev_queue = malloc(sizeof(int) * g->NC);
    ls->in_prev_queue = malloc(sizeof(int) * g->NC);

    ls->fsc = 0;
    ls->clique_status = malloc(sizeof(int) * g->NC);
    ls->fast_set = malloc(sizeof(int) * g->N);

    ls->log_count = 0;
    ls->log = malloc(sizeof(int) * MAX_LOG);
    ls->action = malloc(sizeof(int) * MAX_LOG);

    ls->seed = seed;

    for (int u = 0; u < g->N; u++)
    {
        ls->independent_set[u] = 0;
        ls->fast_set[u] = -1;
    }

    for (int c = 0; c < g->NC; c++)
    {
        ls->queue[c] = c;
        ls->in_queue[c] = 1;
        ls->prev_queue[c] = 0;
        ls->in_prev_queue[c] = 0;

        ls->clique_status[c] = -1;
    }

    return ls;
}

void clique_local_search_free(clique_local_search *ls)
{
    free(ls->independent_set);

    free(ls->queue);
    free(ls->in_queue);
    free(ls->prev_queue);
    free(ls->in_prev_queue);

    free(ls->clique_status);
    free(ls->fast_set);

    free(ls->log);
    free(ls->action);

    free(ls);
}

void clique_local_search_in_order_solution(clique_graph *g, clique_local_search *ls)
{
    for (int u = 0; u < g->N; u++)
    {
        if (ls->independent_set[u])
            continue;

        long long adjacent_weight = 0;
        for (int k = g->V[u]; k < g->V[u + 1]; k++)
            if (ls->clique_status[g->EV[k]] >= 0)
                adjacent_weight += g->W[ls->clique_status[g->EV[k]]];

        if (adjacent_weight < g->W[u])
            clique_local_search_add_vertex(g, ls, u);
    }
}

void clique_local_search_add_vertex(clique_graph *g, clique_local_search *ls, int u)
{
    assert(ls->independent_set[u] == 0);

    ls->action[ls->log_count] = ADD;
    ls->log[ls->log_count] = u;
    ls->log_count++;

    ls->independent_set[u] = 1;
    ls->cost += g->W[u];

    for (int i = g->V[u]; i < g->V[u + 1]; i++)
    {
        int c = g->EV[i];
        long long old_W = 0;
        if (ls->clique_status[c] >= 0)
        {
            old_W = ls->clique_status[c];
            clique_local_search_remove_vertex(g, ls, ls->clique_status[c]);
        }

        ls->clique_status[c] = u;

        if (!ls->in_queue[c] && old_W > g->W[u])
        {
            ls->in_queue[c] = 1;
            ls->queue[ls->queue_count] = c;
            ls->queue_count++;
        }
    }
}

void clique_local_search_undo_add_vertex(clique_graph *g, clique_local_search *ls, int u)
{
    assert(ls->independent_set[u] == 1);

    ls->independent_set[u] = 0;
    ls->cost -= g->W[u];

    for (int i = g->V[u]; i < g->V[u + 1]; i++)
    {
        int c = g->EV[i];
        if (ls->clique_status[c] == u)
            ls->clique_status[c] = -1;
    }
}

void clique_local_search_remove_vertex(clique_graph *g, clique_local_search *ls, int u)
{
    assert(ls->independent_set[u] == 1);

    ls->action[ls->log_count] = REMOVE;
    ls->log[ls->log_count] = u;
    ls->log_count++;

    ls->independent_set[u] = 0;
    ls->cost -= g->W[u];

    for (int i = g->V[u]; i < g->V[u + 1]; i++)
    {
        int c = g->EV[i];
        ls->clique_status[c] = -1;

        if (!ls->in_queue[c])
        {
            ls->in_queue[c] = 1;
            ls->queue[ls->queue_count] = c;
            ls->queue_count++;
        }
    }
}

void clique_local_search_undo_remove_vertex(clique_graph *g, clique_local_search *ls, int u)
{
    assert(ls->independent_set[u] == 0);

    ls->independent_set[u] = 1;
    ls->cost += g->W[u];

    for (int i = g->V[u]; i < g->V[u + 1]; i++)
    {
        int c = g->EV[i];
        ls->clique_status[c] = u;
    }
}

void clique_local_search_shuffle_queue(clique_local_search *ls)
{
    int n = ls->queue_count;
    for (int i = 0; i < n - 1; i++)
    {
        int j = i + rand_r(&ls->seed) / (RAND_MAX / (n - i) + 1);
        int t = ls->queue[j];
        ls->queue[j] = ls->queue[i];
        ls->queue[i] = t;
    }
}

static inline void swap(int **a, int **b)
{
    int *t = *a;
    *a = *b;
    *b = t;
}

void clique_local_search_greedy(clique_graph *g, clique_local_search *ls)
{
    clique_local_search_shuffle_queue(ls);

    int n = ls->queue_count;
    ls->queue_count = 0;
    while (n > 0)
    {
        swap(&ls->queue, &ls->prev_queue);
        swap(&ls->in_queue, &ls->in_prev_queue);

        for (int i = 0; i < n; i++)
        {
            int c = ls->prev_queue[i];
            ls->in_prev_queue[c] = 0;

            for (int j = g->C[c]; j < g->C[c + 1]; j++)
            {
                int u = g->EC[j];

                if (ls->independent_set[u])
                    continue;

                long long adjacent_weight = 0;
                for (int k = g->V[u]; k < g->V[u + 1]; k++)
                    if (ls->clique_status[g->EV[k]] >= 0)
                        adjacent_weight += g->W[ls->clique_status[g->EV[k]]];

                if (adjacent_weight < g->W[u])
                    clique_local_search_add_vertex(g, ls, u);
            }
        }

        clique_local_search_shuffle_queue(ls);

        n = ls->queue_count;
        ls->queue_count = 0;
    }
}

void clique_local_search_explore(clique_graph *g, clique_local_search *ls, int it, int k, int verbose)
{
    int c = 0, t = 0;
    long long best = ls->cost;

    int ref = ls->log_count;

    if (verbose)
    {
        printf("\r%lld %d %d    ", ls->cost, c, t);
        fflush(stdout);
    }

    int q;
    while (c++ < it)
    {
        ls->log_count = ref;
        int p = ls->log_count;

        int u = rand_r(&ls->seed) % g->N;

        if (ls->independent_set[u])
            clique_local_search_remove_vertex(g, ls, u);
        else
            clique_local_search_add_vertex(g, ls, u);

        int r = rand_r(&ls->seed) % k;
        for (int i = 0; i < r && ls->queue_count > 0; i++)
        {
            int c = ls->queue[rand_r(&ls->seed) % ls->queue_count];
            int v = g->EC[g->C[c] + (rand_r(&ls->seed) % (g->C[c + 1] - g->C[c]))];

            if (ls->independent_set[v])
                clique_local_search_remove_vertex(g, ls, v);
            else
                clique_local_search_add_vertex(g, ls, v);
        }

        clique_local_search_greedy(g, ls);

        if (ls->cost > best)
        {
            t++;
            best = ls->cost;
            if (verbose)
            {
                printf("\r%lld %d %d    ", ls->cost, c, t);
                fflush(stdout);
            }
        }
        else if (ls->cost < best)
        {
            clique_local_search_unwind(g, ls, p);
        }
    }
    if (verbose)
        printf("\n");
}

void clique_local_search_unwind(clique_graph *g, clique_local_search *ls, int t)
{
    while (ls->log_count > t)
    {
        ls->log_count--;
        int u = ls->log[ls->log_count];
        int a = ls->action[ls->log_count];

        if (a == ADD)
            clique_local_search_undo_add_vertex(g, ls, u);
        else if (a == REMOVE)
            clique_local_search_undo_remove_vertex(g, ls, u);
    }
}