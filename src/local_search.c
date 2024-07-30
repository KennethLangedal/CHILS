#include "local_search.h"

#include <stdlib.h>

local_search *local_search_init(graph g)
{
    local_search *ls = malloc(sizeof(local_search));

    ls->c = 0;
    ls->IS = malloc(sizeof(int) * g.N);

    ls->qc = g.N;
    ls->a = 1;
    ls->NW = malloc(sizeof(int) * g.N);
    ls->Q = malloc(sizeof(int) * g.N);
    ls->C = malloc(sizeof(int) * g.N);
    ls->_Q = malloc(sizeof(int) * g.N);
    ls->_C = malloc(sizeof(int) * g.N);

    ls->lc = 0;
    ls->A = malloc(sizeof(int) * MAX_LOG);
    ls->L = malloc(sizeof(int) * MAX_LOG);

    for (int u = 0; u < g.N; u++)
    {
        ls->IS[u] = 0;
        ls->NW[u] = 0;
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
    free(ls->Q);
    free(ls->C);
    free(ls->A);
    free(ls->L);
    free(ls);
}

void local_search_add_vertex(graph g, local_search *ls, int u)
{
    if (ls->IS[u])
        return;

    if (ls->a)
    {
        ls->A[ls->lc] = ADD;
        ls->L[ls->lc++] = u;
    }

    ls->IS[u] = 1;
    ls->c += g.W[u];

    for (int i = g.V[u]; i < g.V[u + 1]; i++)
    {
        int v = g.E[i];
        local_search_remove_vertex(g, ls, v);
        ls->NW[v] += g.W[u];

        if (ls->a && !ls->C[v])
        {
            ls->C[v] = 1;
            ls->Q[ls->qc++] = v;
        }
    }
}

void local_search_remove_vertex(graph g, local_search *ls, int u)
{
    if (!ls->IS[u])
        return;

    if (ls->a)
    {
        ls->A[ls->lc] = REMOVE;
        ls->L[ls->lc++] = u;
    }
    ls->IS[u] = 0;
    ls->c -= g.W[u];

    for (int i = g.V[u]; i < g.V[u + 1]; i++)
    {
        int v = g.E[i];
        ls->NW[v] -= g.W[u];

        if (ls->a && !ls->C[v])
        {
            ls->C[v] = 1;
            ls->Q[ls->qc++] = v;
        }
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

            if (ls->NW[u] < g.W[u])
                local_search_add_vertex(g, ls, u);
        }

        local_search_shuffle_queue(ls);

        n = ls->qc;
        ls->qc = 0;
    }
}

void local_search_unwind(graph g, local_search *ls, int t)
{
    ls->a = 0;
    while (ls->lc > t)
    {
        ls->lc--;
        int u = ls->L[ls->lc];
        int a = ls->A[ls->lc];

        if (a == ADD)
            local_search_remove_vertex(g, ls, u);
        else if (a == REMOVE)
            local_search_add_vertex(g, ls, u);
    }
    ls->a = 1;
}

void local_search_shuffle_queue(local_search *ls)
{
    int n = ls->qc;
    for (int i = 0; i < n - 1; i++)
    {
        int j = i + rand() / (RAND_MAX / (n - i) + 1);
        int t = ls->Q[j];
        ls->Q[j] = ls->Q[i];
        ls->Q[i] = t;
    }
}
