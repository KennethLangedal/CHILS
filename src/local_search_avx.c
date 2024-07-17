#include "local_search_avx.h"

#include <stdlib.h>

local_search_avx local_search_init(graph g)
{
    local_search_avx ls;
    ls.N = g.N;
    ls.C = malloc(sizeof(int));
    *ls.C = 0;

    ls.IS = aligned_alloc(32, sizeof(int) * (g.N + 1));
    for (int i = 0; i < g.N + 1; i++)
        ls.IS[i] = 0;

    ls.V = malloc(sizeof(int) * (g.N + 1));
    ls.V[0];
    for (int i = 0; i < g.N; i++)
        ls.V[i + 1] = ls.V[i] + (((g.V[i + 1] - g.V[i]) + 7) & (-8));

    ls.E = aligned_alloc(32, sizeof(int) * ls.V[g.N]);
    ls.W = aligned_alloc(32, sizeof(int) * (g.N + 1));
    for (int i = 0; i < g.N; i++)
    {
        ls.W[i] = g.W[i];
        int j = g.V[i], k = ls.V[i];
        while (j < g.V[i + 1])
            ls.E[k++] = g.E[j++];
        while (k < ls.V[i + 1])
            ls.E[k++] = g.N;
    }
    ls.W[g.N] = 0;

    return ls;
}

void local_search_free(local_search_avx ls)
{
    free(ls.C);
    free(ls.IS);
    free(ls.V);
    free(ls.E);
    free(ls.W);
}

void local_search_add_vertex(local_search_avx ls, int u)
{
    if (ls.IS[u])
        return;

    ls.IS[u] = 0xffffffffu;
    *ls.C += ls.W[u];

    for (int i = ls.V[u]; i < ls.V[u + 1]; i++)
    {
        int v = ls.E[i];
        if (ls.IS[v])
        {
            ls.IS[v] = 0;
            *ls.C -= ls.W[v];
        }
    }
}

void local_search_remove_vertex(local_search_avx ls, int u)
{
    if (!ls.IS[u])
        return;

    ls.IS[u] = 0;
    *ls.C -= ls.W[u];
}

void local_search_greedy(local_search_avx ls, int *order)
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