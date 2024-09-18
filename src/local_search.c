#include "local_search.h"

#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <omp.h>

#define REMOVE 0
#define ADD 1

#define MAX_LOG 10000000
#define MAX_GUESS 100
#define MAX_TWO_ONE_DEGREE 256

local_search *local_search_init(graph *g, unsigned int seed)
{
    local_search *ls = malloc(sizeof(local_search));

    ls->cost = 0;
    ls->independent_set = malloc(sizeof(int) * g->N);

    ls->queue_count = g->N;
    ls->queue = malloc(sizeof(int) * g->N);
    ls->in_queue = malloc(sizeof(int) * g->N);
    ls->prev_queue = malloc(sizeof(int) * g->N);
    ls->in_prev_queue = malloc(sizeof(int) * g->N);

    ls->adjacent_weight = malloc(sizeof(long long) * g->N);
    ls->tabu = malloc(sizeof(int) * g->N);
    ls->tightness = malloc(sizeof(int) * g->N);
    ls->temp = malloc(sizeof(int) * g->N * 2);
    ls->mask = malloc(sizeof(int) * g->N);

    ls->log_count = 0;
    ls->log = malloc(sizeof(int) * MAX_LOG);
    ls->action = malloc(sizeof(int) * MAX_LOG);

    ls->seed = seed;

    for (int u = 0; u < g->N; u++)
    {
        ls->independent_set[u] = 0;
        ls->queue[u] = u;
        ls->in_queue[u] = 1;
        ls->prev_queue[u] = 0;
        ls->in_prev_queue[u] = 0;

        ls->adjacent_weight[u] = 0;
        ls->tabu[u] = 0;
        ls->tightness[u] = 0;
        ls->temp[u] = 0;
        ls->mask[u] = 0;
    }

    return ls;
}

void local_search_free(local_search *ls)
{
    free(ls->independent_set);

    free(ls->queue);
    free(ls->in_queue);
    free(ls->prev_queue);
    free(ls->in_prev_queue);

    free(ls->adjacent_weight);
    free(ls->tabu);
    free(ls->tightness);
    free(ls->temp);
    free(ls->mask);

    free(ls->log);
    free(ls->action);

    free(ls);
}

void local_search_in_order_solution(graph *g, local_search *ls)
{
    for (int u = 0; u < g->N; u++)
    {
        if (!ls->tabu[u] && ls->adjacent_weight[u] < g->W[u])
            local_search_add_vertex(g, ls, u);
    }
}

void local_search_add_vertex(graph *g, local_search *ls, int u)
{
    assert(!ls->independent_set[u] && !ls->tabu[u]);

    ls->action[ls->log_count] = ADD;
    ls->log[ls->log_count] = u;
    ls->log_count++;

    ls->independent_set[u] = 1;
    ls->cost += g->W[u];

    for (int i = g->V[u]; i < g->V[u + 1]; i++)
    {
        int v = g->E[i];
        if (ls->independent_set[v])
            local_search_remove_vertex(g, ls, v);

        ls->adjacent_weight[v] += g->W[u];
        ls->tightness[v]++;
    }
}

void local_search_undo_add_vertex(graph *g, local_search *ls, int u)
{
    assert(ls->independent_set[u] && !ls->tabu[u]);

    ls->independent_set[u] = 0;
    ls->cost -= g->W[u];

    for (int i = g->V[u]; i < g->V[u + 1]; i++)
    {
        int v = g->E[i];
        ls->adjacent_weight[v] -= g->W[u];
        ls->tightness[v]--;
    }
}

void local_search_remove_vertex(graph *g, local_search *ls, int u)
{
    assert(ls->independent_set[u] && !ls->tabu[u]);

    ls->action[ls->log_count] = REMOVE;
    ls->log[ls->log_count] = u;
    ls->log_count++;

    ls->independent_set[u] = 0;
    ls->cost -= g->W[u];

    if (!ls->in_queue[u])
    {
        ls->in_queue[u] = 1;
        ls->queue[ls->queue_count] = u;
        ls->queue_count++;
    }

    for (int i = g->V[u]; i < g->V[u + 1]; i++)
    {
        int v = g->E[i];
        ls->adjacent_weight[v] -= g->W[u];
        ls->tightness[v]--;

        if (!ls->in_queue[v] && !ls->tabu[v])
        {
            ls->in_queue[v] = 1;
            ls->queue[ls->queue_count] = v;
            ls->queue_count++;
        }
    }
}

void local_search_undo_remove_vertex(graph *g, local_search *ls, int u)
{
    assert(!ls->independent_set[u] && !ls->tabu[u]);

    ls->independent_set[u] = 1;
    ls->cost += g->W[u];

    for (int i = g->V[u]; i < g->V[u + 1]; i++)
    {
        int v = g->E[i];
        ls->adjacent_weight[v] += g->W[u];
        ls->tightness[v]++;
    }
}

void local_search_lock_vertex(graph *g, local_search *ls, int u)
{
    ls->tabu[u]++;
    if (ls->independent_set[u])
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
            ls->tabu[g->E[i]]++;
}

void local_search_unlock_vertex(graph *g, local_search *ls, int u)
{
    ls->tabu[u]--;
    if (!ls->in_queue[u] && !ls->tabu[u])
    {
        ls->in_queue[u] = 1;
        ls->queue[ls->queue_count] = u;
        ls->queue_count++;
    }
    if (!ls->independent_set[u])
        return;

    for (int i = g->V[u]; i < g->V[u + 1]; i++)
    {
        int v = g->E[i];
        ls->tabu[v]--;

        if (!ls->in_queue[v] && !ls->tabu[v])
        {
            ls->in_queue[v] = 1;
            ls->queue[ls->queue_count] = v;
            ls->queue_count++;
        }
    }
}

void local_search_two_one(graph *g, local_search *ls, int u)
{
    assert(ls->independent_set[u] && !ls->tabu[u]);

    int adjacent_count = 0;
    for (int i = g->V[u]; i < g->V[u + 1]; i++)
        if (ls->tightness[g->E[i]] == 1 && !ls->tabu[g->E[i]])
            ls->temp[adjacent_count++] = g->E[i];

    if (adjacent_count < 2)
        return;

    for (int i = 0; i < adjacent_count; i++)
    {
        int v = ls->temp[i];

        int i1 = 0, i2 = g->V[v];
        while (i1 < adjacent_count && i2 < g->V[v + 1])
        {
            int w1 = ls->temp[i1], w2 = g->E[i2];
            if (w1 > w2)
                i2++;
            else if (w1 == w2)
                i1++, i2++;
            else if (w1 == v)
                i1++;
            else if (g->W[w1] + g->W[v] > g->W[u])
            {
                local_search_add_vertex(g, ls, v);
                local_search_add_vertex(g, ls, w1);
                return;
            }
            else
                i2++;
        }
    }
}

void local_search_aap(graph *g, local_search *ls, int u)
{
    assert(ls->independent_set[u] || ls->tightness[u] == 1);

    int current = -1, candidate_size = 0;
    ls->temp[candidate_size++] = u;

    if (ls->independent_set[u])
    {
        current = u;
    }
    else
    {
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (ls->independent_set[v])
            {
                current = v;
                ls->temp[candidate_size++] = v;
                break;
            }
        }

        if (ls->tabu[current])
            return;
        ls->mask[u] = 1;
    }

    int found = 1;
    while (found)
    {
        found = 0;
        long long best = INT_MIN;
        int to_add, to_remove;

        for (int i = g->V[current]; i < g->V[current + 1]; i++)
        {
            int v = g->E[i];
            if (ls->tightness[v] != 2 || ls->mask[v] || ls->tabu[v])
                continue;

            int valid = 1, next = -1;
            for (int j = g->V[v]; j < g->V[v + 1] && valid; j++)
            {
                int w = g->E[j];
                if (w == current)
                    continue;

                if (ls->mask[w])
                    valid = 0;
                else if (ls->independent_set[w])
                    next = w;
            }

            long long gain = (rand_r(&ls->seed) % (1 << 10)) - (1 << 9);
            if (valid && !ls->tabu[next] && (g->W[v] - g->W[next]) + gain > best)
            {
                to_add = v;
                to_remove = next;
                best = (g->W[v] - g->W[next]) + gain;
                found = 1;
            }
        }

        if (found)
        {
            ls->temp[candidate_size++] = to_add;
            ls->mask[to_add] = 1;
            ls->temp[candidate_size++] = to_remove;
            current = to_remove;
        }
    }

    long long diff = 0, best = LLONG_MIN;
    int best_position = 0;
    int to_add = 0;

    for (int i = 0; i < candidate_size; i++)
    {
        int v = ls->temp[i];
        if (ls->independent_set[v] && !ls->mask[v])
        {
            diff -= g->W[v];
            ls->mask[v] = 1;
        }
        else if (!ls->independent_set[v])
        {
            diff += g->W[v];
            ls->temp[g->N + to_add++] = v;
        }

        if (ls->independent_set[v] && diff > best)
        {
            best = diff;
            best_position = to_add;
        }
    }

    if (best <= 0)
        best_position = to_add;

    for (int i = 0; i < best_position; i++)
    {
        int v = ls->temp[g->N + i];
        local_search_add_vertex(g, ls, v);
    }

    for (int i = 0; i < candidate_size; i++)
        ls->mask[ls->temp[i]] = 0;
}

void local_search_shuffle_queue(local_search *ls)
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

void swap(int **a, int **b)
{
    int *t = *a;
    *a = *b;
    *b = t;
}

void local_search_greedy(graph *g, local_search *ls)
{
    local_search_shuffle_queue(ls);

    int n = ls->queue_count;
    ls->queue_count = 0;
    while (n > 0)
    {
        swap(&ls->queue, &ls->prev_queue);
        swap(&ls->in_queue, &ls->in_prev_queue);

        for (int i = 0; i < n; i++)
        {
            int u = ls->prev_queue[i];
            ls->in_prev_queue[u] = 0;

            if (ls->tabu[u])
                continue;

            if (!ls->independent_set[u] && ls->adjacent_weight[u] < g->W[u])
                local_search_add_vertex(g, ls, u);
            else if (ls->independent_set[u] && g->V[u + 1] - g->V[u] < MAX_TWO_ONE_DEGREE)
                local_search_two_one(g, ls, u);
        }

        local_search_shuffle_queue(ls);

        n = ls->queue_count;
        ls->queue_count = 0;
    }
}

void local_search_explore(graph *g, local_search *ls, double tl, int verbose, long long offset)
{
    int c = 0, q = 0;

    long long best = ls->cost;

    if (verbose)
    {
        printf("Running local search for %.2lf seconds\n", tl);
        printf("\r%lld %.2lf  ", ls->cost + offset, 0.0);
        fflush(stdout);
    }

    double start, end;
    start = omp_get_wtime();

    if (ls->cost == 0)
        local_search_in_order_solution(g, ls);

    local_search_greedy(g, ls);

    while (1)
    {
        if ((c++ & ((1 << 10) - 1)) == 0)
        {
            c = 0;
            end = omp_get_wtime();
            double elapsed = end - start;
            if (elapsed > tl)
                break;
        }

        ls->log_count = 0;

        int u = rand_r(&ls->seed) % g->N;
        q = 0;
        while (q++ < MAX_GUESS && ls->tabu[u])
            u = rand_r(&ls->seed) % g->N;

        if (ls->tabu[u])
            continue;

        if (ls->independent_set[u] || ls->tightness[u] == 1)
        {
            local_search_aap(g, ls, u);

            local_search_greedy(g, ls);
        }
        else
        {
            local_search_add_vertex(g, ls, u);
            local_search_lock_vertex(g, ls, u);

            // int to_remove = __builtin_clz(((unsigned int)rand_r(&ls->seed)) + 1);
            int to_remove = rand_r(&ls->seed) & 15;
            for (int i = 0; i < to_remove && ls->queue_count > 0 && ls->queue_count < (1 << 12); i++)
            {
                int v = ls->queue[rand_r(&ls->seed) % ls->queue_count];
                q = 0;
                while (q++ < MAX_GUESS && ls->tabu[v])
                    v = ls->queue[rand_r(&ls->seed) % ls->queue_count];

                if (ls->tabu[v])
                    continue;

                if (ls->independent_set[v])
                    local_search_remove_vertex(g, ls, v);
                else
                    local_search_add_vertex(g, ls, v);
            }

            local_search_greedy(g, ls);
            local_search_unlock_vertex(g, ls, u);
            local_search_greedy(g, ls);
        }

        if (ls->cost > best)
        {
            best = ls->cost;
            ls->log_count = 0;
            if (verbose)
            {
                end = omp_get_wtime();
                double elapsed = end - start;

                printf("\r%lld %.2lf  ", ls->cost + offset, elapsed);
                fflush(stdout);
            }
        }
        if (ls->cost < best)
        {
            local_search_unwind(g, ls, 0);
        }
    }
    if (verbose)
        printf("\n");
}

void local_search_unwind(graph *g, local_search *ls, int t)
{
    while (ls->log_count > t)
    {
        ls->log_count--;
        int u = ls->log[ls->log_count];
        int a = ls->action[ls->log_count];

        if (a == ADD)
            local_search_undo_add_vertex(g, ls, u);
        else if (a == REMOVE)
            local_search_undo_remove_vertex(g, ls, u);
    }
}
