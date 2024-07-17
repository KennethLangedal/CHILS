#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "local_search.h"
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

    local_search ls = local_search_init(g);

    int *order = malloc(sizeof(int) * g.N);
    for (int i = 0; i < g.N; i++)
        order[i] = i;

    shuffle(order, g.N);

    local_search_greedy(g, ls, order);
    // local_search_k_one(g, ls, order);
    // local_search_k_c(g, ls, order);

    int C = *ls.C;
    int *IS = malloc(sizeof(int) * g.N);
    for (int i = 0; i < g.N; i++)
        IS[i] = ls.IS[i];

    int k = 0;
    while (k++ < 1000)
    {
        for (int i = 0; i < 1; i++)
        {
            int u = rand() % g.N;
            if (ls.IS[u])
                local_search_remove_vertex(g, ls, u);
            else
                local_search_add_vertex(g, ls, u);
        }

        shuffle(order, g.N);

        local_search_greedy(g, ls, order);
        // local_search_k_one(g, ls, order);
        // local_search_k_c(g, ls, order);

        if (*ls.C > C)
        {
            C = *ls.C;
            for (int i = 0; i < g.N; i++)
                IS[i] = ls.IS[i];

            valid = mwis_validate(g, ls.IS);
            printf("%d %d\n", *ls.C, valid);
        }
        else if (*ls.C < C)
        {
            *ls.C = C;
            for (int i = 0; i < g.N; i++)
                ls.IS[i] = IS[i];
        }
    }

    graph_free(g);
    local_search_free(ls);
    free(order);

    return 0;
}