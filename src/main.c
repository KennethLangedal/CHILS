#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "local_search.h"
#include "local_search_avx.h"
#include "mwis.h"

void shuffle(int *array, int n)
{
    for (int i = 0; i < n - 1; i++)
    {
        int j = i + rand() / (RAND_MAX / (n - i) + 1);
        int t = array[j];
        array[j] = array[i];
        array[i] = t;
    }
}

int main(int argc, char **argv)
{
    graph g = graph_parse(stdin);

    int valid = graph_validate(g.N, g.V, g.E);
    printf("Graph valid %d |V|=%d, |E|=%d\n", valid, g.N, g.V[g.N] / 2);

    local_search ls = local_search_init(g), best = local_search_init(g);
    // local_search_avx ls = local_search_avx_init(g);

    int *order = malloc(sizeof(int) * g.N);
    for (int i = 0; i < g.N; i++)
        order[i] = i;

    int *order_pos = malloc(sizeof(int) * g.N);
    for (int i = 0; i < g.N; i++)
        order_pos[i] = i;

    int *in_count = malloc(sizeof(int) * g.N);
    for (int i = 0; i < g.N; i++)
        in_count[i] = 0;

    // shuffle(order, g.N);

    local_search_greedy(g, ls, order);
    // local_search_k_one(g, ls, order);
    // local_search_k_c(g, ls, order);

    local_search_copy(g, ls, best);

    int k = 0, t = 0;
    while (k++ < 10000000)
    {
        for (int i = 0; i < (rand() % 8) + 1; i++)
        {
            int u = rand() % g.N;
            if (ls.IS[u])
                local_search_remove_vertex(g, ls, u);
            else
                local_search_add_vertex(g, ls, u);

            int p = order_pos[u];
            int v = order[g.N - 1 - i];
            order[g.N - 1 - i] = u;
            order_pos[u] = g.N - 1 - i;
            order[p] = v;
            order_pos[v] = p;
        }

        // if ((k % 1000) == 0)
        //     shuffle(order, g.N);

        local_search_greedy(g, ls, order);

        // ls.tabu[u] = 0;
        // ls.tabu[v] = 0;

        if (*ls.C > *best.C)
        {
            local_search_copy(g, ls, best);

            t++;
            // for (int i = 0; i < g.N; i++)
            //     if (ls.IS[i])
            //         in_count[i]++;

            // valid = mwis_validate(g, ls.IS);
            printf("%d %d %d\n", *ls.C, k, t);
        }
        else if (*ls.C < *best.C)
        {
            local_search_copy(g, best, ls);
        }
    }

    // int rc = 0;
    // for (int i = 0; i < g.N; i++)
    // {
    //     if (in_count[i] == 0)
    //         local_search_lock_out_vertex(g, ls, i);
    //     else if (in_count[i] == t)
    //         local_search_lock_in_vertex(g, ls, i);
    //     else
    //         continue;

    //     rc++;
    // }

    // printf("Removed %d vertices\n", rc);

    // k = 0, t = 0;
    // while (k++ < 10000)
    // {
    //     for (int i = 0; i < 1; i++)
    //     {
    //         int u = rand() % g.N;
    //         while (ls.tabu[u])
    //             u = rand() % g.N;

    //         if (ls.IS[u])
    //             local_search_remove_vertex(g, ls, u);
    //         else
    //             local_search_add_vertex(g, ls, u);
    //     }

    //     shuffle(order, g.N);

    //     local_search_greedy(g, ls, order);
    //     // local_search_k_one(g, ls, order);
    //     // local_search_k_c(g, ls, order);

    //     if (*ls.C >= C)
    //     {
    //         C = *ls.C;
    //         for (int i = 0; i < g.N; i++)
    //             IS[i] = ls.IS[i];

    //         t++;
    //         for (int i = 0; i < g.N; i++)
    //             if (IS[i])
    //                 in_count[i]++;

    //         valid = mwis_validate(g, ls.IS);
    //         printf("%d %d\n", *ls.C, valid);
    //     }
    //     else if (*ls.C < C)
    //     {
    //         *ls.C = C;
    //         for (int i = 0; i < g.N; i++)
    //             ls.IS[i] = IS[i];
    //     }
    // }

    graph_free(g);
    local_search_free(ls);
    local_search_free(best);
    free(order);

    return 0;
}