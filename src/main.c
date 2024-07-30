#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "local_search.h"
#include "mwis.h"

int *run_ls(graph g, int rp, int k, int offset)
{
    local_search *ls = local_search_init(g);

    local_search_greedy(g, ls);

    int best = ls->c;
    int c = 0, t = 0;

    while (c++ < k)
    {
        for (int j = 0; j < g.N; j++)
        {
            int p = ls->lc;

            int u = j; // rand() % g.N;
            if (ls->IS[u])
                local_search_remove_vertex(g, ls, u);
            else
                local_search_add_vertex(g, ls, u);

            int r = rand() % rp;
            for (int i = 0; i < r && ls->qc > 0; i++)
            {
                u = ls->Q[rand() % ls->qc];

                if (ls->IS[u])
                    local_search_remove_vertex(g, ls, u);
                else
                    local_search_add_vertex(g, ls, u);
            }

            local_search_greedy(g, ls);

            if (ls->c > best)
            {
                t++;
                best = ls->c;
                printf("\r%d %d %d %d", ls->c + offset, c, t, ls->lc);
                fflush(stdout);
            }
            else if (ls->c <= best)
            {
                local_search_unwind(g, ls, p);
            }
        }
    }
    printf("\n");

    int *S = malloc(sizeof(int) * g.N);
    for (int i = 0; i < g.N; i++)
        S[i] = ls->IS[i];

    local_search_free(ls);

    return S;
}

int main(int argc, char **argv)
{
    graph g = graph_parse(stdin);

    int valid = graph_validate(g.N, g.V, g.E);
    printf("Graph valid %d |V|=%d, |E|=%d\n", valid, g.N, g.V[g.N] / 2);

    int *C = malloc(sizeof(int) * g.N);
    for (int i = 0; i < g.N; i++)
        C[i] = 0;

    for (int i = 0; i < 10; i++)
    {
        int *S = run_ls(g, 64, 10000000, 0);
        int in = 0, out = 0;
        for (int j = 0; j < g.N; j++)
        {
            C[j] += S[j];
            if (C[j] == 0)
                out++;
            else if (C[j] == i + 1)
                in++;
        }
        printf("%d %d / %d\n", in, out, g.N);
        free(S);
    }

    int *mask = malloc(sizeof(int) * g.N);
    int offset = 0;
    for (int i = 0; i < g.N; i++)
        mask[i] = 1;

    for (int u = 0; u < g.N; u++)
    {
        if (C[u] < 10)
            continue;
        offset += g.W[u];
        mask[u] = 0;
        for (int i = g.V[u]; i < g.V[u + 1]; i++)
        {
            int v = g.E[i];
            mask[v] = 0;
        }
    }

    graph sg = graph_subgraph(g, mask, C);

    for (int i = 0; i < 10; i++)
    {
        int *S = run_ls(sg, 32, 1000000, offset);
        free(S);
    }

    graph_free(g);

    return 0;
}