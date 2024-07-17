#include "local_search.h"

#include <stdlib.h>

local_search local_search_init(graph g)
{
    local_search ls;
    ls.C = malloc(sizeof(int));
    ls.IS = malloc(sizeof(int) * g.N);

    *ls.C = 0;

    for (int u = 0; u < g.N; u++)
        ls.IS[u] = 0;

    return ls;
}

void local_search_free(local_search ls)
{
    free(ls.C);
    free(ls.IS);
}

void local_search_add_vertex(graph g, local_search ls, int u)
{
    if (ls.IS[u])
        return;

    ls.IS[u] = 1;
    *ls.C += g.W[u];

    for (int i = g.V[u]; i < g.V[u + 1]; i++)
    {
        int v = g.E[i];
        if (ls.IS[v])
        {
            ls.IS[v] = 0;
            *ls.C -= g.W[v];
        }
    }
}

void local_search_remove_vertex(graph g, local_search ls, int u)
{
    if (!ls.IS[u])
        return;

    ls.IS[u] = 0;
    *ls.C -= g.W[u];
}

void local_search_greedy(graph g, local_search ls, int *order)
{
    int imp = 1;
    while (imp)
    {
        imp = 0;
        for (int i = 0; i < g.N; i++)
        {
            int u = order[i];
            if (ls.IS[u])
                continue;

            int nw = 0;
            for (int j = g.V[u]; j < g.V[u + 1]; j++)
            {
                int v = g.E[j];
                if (!ls.IS[v])
                    continue;

                nw += g.W[v];
                if (nw >= g.W[u])
                    break;
            }

            if (nw < g.W[u])
            {
                local_search_add_vertex(g, ls, u);
                imp = 1;
            }
        }
    }
}

void local_search_k_one(graph g, local_search ls, int *order)
{
    int imp = 1;

    graph _g;

    _g.N = 0;
    _g.V = malloc(sizeof(int) * (g.N + 1));
    _g.E = malloc(sizeof(int) * g.V[g.N]);
    _g.W = malloc(sizeof(int) * g.N);

    int *R = malloc(sizeof(int) * g.N);
    int *tmp = malloc(sizeof(int) * g.N);
    int *mask = malloc(sizeof(int) * g.N);

    while (imp)
    {
        for (int i = 0; i < g.N; i++)
            mask[i] = -1;

        imp = 0;
        for (int i = 0; i < g.N; i++)
        {
            int u = order[i];
            if (!ls.IS[u])
                continue;

            int nc = 0, nw = 0;
            for (int j = g.V[u]; j < g.V[u + 1]; j++)
            {
                int v = g.E[j];
                int val = 1;
                for (int k = g.V[v]; k < g.V[v + 1]; k++)
                {
                    int w = g.E[k];
                    if (w != u && ls.IS[w])
                    {
                        val = 0;
                        break;
                    }
                }

                if (!val)
                    continue;

                R[v] = nc;
                tmp[nc] = v;
                mask[v] = u;
                nc++;
                nw += g.W[v];
            }

            if (nc < 2 || nw < g.W[u])
                continue;

            // Construct subgraph
            _g.N = nc;
            int ei = 0;
            for (int j = 0; j < _g.N; j++)
            {
                int v = tmp[j];
                _g.V[j] = ei;
                _g.W[j] = g.W[v];

                for (int k = g.V[v]; k < g.V[v + 1]; k++)
                {
                    int w = g.E[k];
                    if (mask[w] == u)
                        _g.E[ei++] = R[w];
                }
            }
            _g.V[_g.N] = ei;

            for (int j = 0; j < _g.N; j++)
                R[j] = j;

            local_search _ls = local_search_init(_g);
            local_search_greedy(_g, _ls, R);

            // printf("%d vs %d\n", g.W[u], *_ls.C);
            if (g.W[u] < *_ls.C)
            {
                imp = 1;
                local_search_remove_vertex(g, ls, u);
                for (int j = 0; j < nc; j++)
                {
                    int v = tmp[j];
                    if (_ls.IS[j])
                        local_search_add_vertex(g, ls, v);
                }
            }

            local_search_free(_ls);
        }
    }

    free(R);
    free(tmp);
    free(mask);

    graph_free(_g);
}

void local_search_k_c(graph g, local_search ls, int *order)
{
    int imp = 1;

    graph _g;

    _g.N = 0;
    _g.V = malloc(sizeof(int) * (g.N + 1));
    _g.E = malloc(sizeof(int) * g.V[g.N]);
    _g.W = malloc(sizeof(int) * g.N);

    int *R = malloc(sizeof(int) * g.N);
    int *tmp = malloc(sizeof(int) * g.N);
    int *mask = malloc(sizeof(int) * g.N);

    while (imp)
    {
        for (int i = 0; i < g.N; i++)
            mask[i] = -1;

        imp = 0;
        for (int i = 0; i < g.N; i++)
        {
            int u = order[i];
            if (!ls.IS[u])
                continue;

            int nc = 0, nw = 0, rw = g.W[u];
            for (int j = g.V[u]; j < g.V[u + 1]; j++)
            {
                int v = g.E[j];
                for (int k = g.V[v]; k < g.V[v + 1]; k++)
                {
                    int w = g.E[k];
                    if (w != u && ls.IS[w] && mask[w] != u)
                    {
                        R[w] = nc;
                        tmp[nc] = w;
                        mask[w] = u;
                        nc++;
                        rw += g.W[w];
                    }
                }

                R[v] = nc;
                tmp[nc] = v;
                mask[v] = u;
                nc++;
                nw += g.W[v];
            }

            if (nc < 1 || nw < g.W[u])
                continue;

            // Construct subgraph
            _g.N = nc;
            int ei = 0;
            for (int j = 0; j < _g.N; j++)
            {
                int v = tmp[j];
                _g.V[j] = ei;
                _g.W[j] = g.W[v];

                for (int k = g.V[v]; k < g.V[v + 1]; k++)
                {
                    int w = g.E[k];
                    if (mask[w] == u)
                        _g.E[ei++] = R[w];
                }
            }
            _g.V[_g.N] = ei;

            for (int j = 0; j < _g.N; j++)
                R[j] = j;

            local_search _ls = local_search_init(_g);
            local_search_greedy(_g, _ls, R);

            if (rw < *_ls.C)
            {
                imp = 1;
                for (int j = 0; j < nc; j++)
                {
                    int v = tmp[j];
                    if (_ls.IS[j])
                        local_search_add_vertex(g, ls, v);
                }
                // printf("%d vs %d (%d)\n", rw, *_ls.C, *ls.C);
            }

            local_search_free(_ls);
        }
    }

    free(R);
    free(tmp);
    free(mask);

    graph_free(_g);
}