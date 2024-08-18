#include "local_search.h"

#include <stdlib.h>
#include <assert.h>

local_search *local_search_init(graph g, unsigned int seed)
{
    local_search *ls = malloc(sizeof(local_search));

    ls->c = 0;
    ls->IS = malloc(sizeof(int) * g.N);

    ls->qc = g.N;
    ls->NW = malloc(sizeof(int) * g.N);
    ls->T = malloc(sizeof(int) * g.N);
    ls->_T = malloc(sizeof(int) * g.N);
    ls->P = malloc(sizeof(int) * g.N);
    ls->_P = malloc(sizeof(int) * g.N);
    ls->Q = malloc(sizeof(int) * g.N);
    ls->C = malloc(sizeof(int) * g.N);
    ls->_Q = malloc(sizeof(int) * g.N);
    ls->_C = malloc(sizeof(int) * g.N);

    ls->lc = 0;
    ls->A = malloc(sizeof(int) * MAX_LOG);
    ls->L = malloc(sizeof(int) * MAX_LOG);

    ls->seed = seed;

    for (int u = 0; u < g.N; u++)
    {
        ls->IS[u] = 0;
        ls->NW[u] = 0;
        ls->T[u] = 0;
        ls->_T[u] = 0;
        ls->P[u] = 0;
        ls->_P[u] = 0;
        ls->Q[u] = u;
        ls->C[u] = 1;
        ls->_Q[u] = 0;
        ls->_C[u] = 0;
    }

    return ls;
}

void local_search_free(local_search *ls)
{
    free(ls->IS);
    free(ls->NW);
    free(ls->T);
    free(ls->_T);
    free(ls->P);
    free(ls->_P);
    free(ls->Q);
    free(ls->C);
    free(ls->_Q);
    free(ls->_C);
    free(ls->A);
    free(ls->L);
    free(ls);
}

void local_search_add_vertex(graph g, local_search *ls, int u)
{
    assert(ls->IS[u] == 0);

    ls->A[ls->lc] = ADD;
    ls->L[ls->lc++] = u;

    ls->IS[u] = 1;
    ls->c += g.W[u];

    if (!ls->C[u] && !ls->T[u])
    {
        ls->C[u] = 1;
        ls->Q[ls->qc++] = u;
    }

    for (int i = g.V[u]; i < g.V[u + 1]; i++)
    {
        int v = g.E[i];
        if (ls->IS[v])
            local_search_remove_vertex(g, ls, v);
        ls->NW[v] += g.W[u];
        ls->P[v]++;
    }
}

void local_search_undo_add_vertex(graph g, local_search *ls, int u)
{
    assert(ls->IS[u] == 1);

    ls->IS[u] = 0;
    ls->c -= g.W[u];

    for (int i = g.V[u]; i < g.V[u + 1]; i++)
    {
        int v = g.E[i];
        ls->NW[v] -= g.W[u];
        ls->P[v]--;
    }
}

void local_search_remove_vertex(graph g, local_search *ls, int u)
{
    assert(ls->IS[u] == 1);

    ls->A[ls->lc] = REMOVE;
    ls->L[ls->lc++] = u;

    ls->IS[u] = 0;
    ls->c -= g.W[u];

    if (!ls->C[u] && !ls->T[u])
    {
        ls->C[u] = 1;
        ls->Q[ls->qc++] = u;
    }

    for (int i = g.V[u]; i < g.V[u + 1]; i++)
    {
        int v = g.E[i];
        ls->NW[v] -= g.W[u];
        ls->P[v]--;

        if (!ls->C[v] && !ls->T[v])
        {
            ls->C[v] = 1;
            ls->Q[ls->qc++] = v;
        }
    }
}

void local_search_undo_remove_vertex(graph g, local_search *ls, int u)
{
    assert(ls->IS[u] == 0);

    ls->IS[u] = 1;
    ls->c += g.W[u];

    for (int i = g.V[u]; i < g.V[u + 1]; i++)
    {
        int v = g.E[i];
        ls->NW[v] += g.W[u];
        ls->P[v]++;
    }
}

void local_search_lock_vertex(graph g, local_search *ls, int u)
{
    ls->T[u] = 1;
    if (!ls->IS[u])
        return;

    for (int i = g.V[u]; i < g.V[u + 1]; i++)
        ls->T[u] = 1;
}

void local_search_unlock_vertex(graph g, local_search *ls, int u)
{
    ls->T[u] = 0;
    if (!ls->IS[u])
        return;

    for (int i = g.V[u]; i < g.V[u + 1]; i++)
        ls->T[u] = 0;
}

void local_search_two_one(graph g, local_search *ls, int u)
{
    if (!ls->IS[u])
        return;

    int nc = 0;
    for (int i = g.V[u]; i < g.V[u + 1]; i++)
        if (ls->NW[g.E[i]] == g.W[u])
            ls->_T[nc++] = g.E[i];

    if (nc < 2)
        return;

    for (int i = 0; i < nc; i++)
    {
        int v = ls->_T[i];

        int i1 = 0, i2 = g.V[v];
        while (i1 < nc && i2 < g.V[v + 1])
        {
            int w1 = ls->_T[i1], w2 = g.E[i2];
            if (w1 > w2)
                i2++;
            else if (w1 == w2)
                i1++, i2++;
            else if (w1 == v)
                i1++;
            else if (g.W[w1] + g.W[v] > g.W[u])
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

void local_search_aap(graph g, local_search *ls, int u)
{
    assert(ls->IS[u] == 0 && ls->P[u] == 1);

    int nc = 0;
    ls->_T[nc++] = u;
    ls->_P[u] = 1;

    int c = -1;

    for (int i = g.V[u]; i < g.V[u + 1]; i++)
    {
        int v = g.E[i];
        if (ls->IS[v])
        {
            c = v;
            ls->_T[nc++] = v;
            break;
        }
    }

    int found = 1;
    while (found)
    {
        found = 0;
        int best = -9999999;
        int w, x;

        for (int i = g.V[c]; i < g.V[c + 1]; i++)
        {
            int v = g.E[i];
            if (ls->P[v] == 2 && !ls->_P[v])
            {
                int val = 1, next = -1;
                for (int j = g.V[v]; j < g.V[v + 1]; j++)
                {
                    int w = g.E[j];
                    if (w == c)
                        continue;

                    if (ls->_P[w])
                    {
                        val = 0;
                        break;
                    }
                    if (ls->IS[w])
                        next = w;
                }

                int gain = (rand_r(&ls->seed) % (1 << 16)) - (1 << 15);
                if (val && (g.W[v] - g.W[next]) + gain > best)
                {
                    w = v;
                    x = next;
                    best = (g.W[v] - g.W[next]) + gain;
                    found = 1;
                }
            }
        }

        if (found)
        {
            ls->_T[nc++] = w;
            ls->_P[w] = 1;
            ls->_T[nc++] = x;
            c = x;
        }
    }

    long long diff = 0, best = 0;
    int best_n = 0;
    int to_add = 0;

    for (int i = 0; i < nc; i++)
    {
        int v = ls->_T[i];
        if (ls->IS[v] && !ls->_P[v])
        {
            diff -= g.W[v];
            ls->_P[v] = 1;
        }
        else if (!ls->IS[v])
        {
            diff += g.W[v];
            ls->_T[to_add++] = v;
        }

        if (ls->IS[v] && diff > best)
        {
            best = diff;
            best_n = to_add;
        }
    }

    for (int i = 0; i < best_n; i++)
    {
        int v = ls->_T[i];
        local_search_add_vertex(g, ls, v);
    }

    for (int i = 0; i < nc; i++)
        ls->_P[ls->_T[i]] = 0;

    return;

    if (best > 0)
        printf("\n%d %lld\n", best_n, best);
}

void local_search_shuffle_queue(local_search *ls)
{
    int n = ls->qc;
    for (int i = 0; i < n - 1; i++)
    {
        int j = i + rand_r(&ls->seed) / (RAND_MAX / (n - i) + 1);
        int t = ls->Q[j];
        ls->Q[j] = ls->Q[i];
        ls->Q[i] = t;
    }
}

void swap(int **a, int **b)
{
    int *t = *a;
    *a = *b;
    *b = t;
}

void local_search_greedy(graph g, local_search *ls)
{
    local_search_shuffle_queue(ls);

    int n = ls->qc;
    ls->qc = 0;
    while (n > 0)
    {
        swap(&ls->Q, &ls->_Q);
        swap(&ls->C, &ls->_C);

        for (int i = 0; i < n; i++)
        {
            int u = ls->_Q[i];
            ls->_C[u] = 0;

            if (ls->T[u])
                continue;

            if (!ls->IS[u] && ls->NW[u] < g.W[u])
                local_search_add_vertex(g, ls, u);
            else if (ls->IS[u])
                local_search_two_one(g, ls, u);
        }

        local_search_shuffle_queue(ls);

        n = ls->qc;
        ls->qc = 0;
    }
}

void local_search_unwind(graph g, local_search *ls, int t)
{
    while (ls->lc > t)
    {
        ls->lc--;
        int u = ls->L[ls->lc];
        int a = ls->A[ls->lc];

        if (a == ADD)
            local_search_undo_add_vertex(g, ls, u);
        else if (a == REMOVE)
            local_search_undo_remove_vertex(g, ls, u);
    }
}
