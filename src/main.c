#include <stdio.h>
#include <stdlib.h>
// #include <omp.h>

#include "graph.h"
#include "local_search.h"
#include "reductions.h"
#include "kernelization.h"
#include "mwis.h"

// static inline int compare(const void *a, const void *b)
// {
//     return *(int *)a - *(int *)b;
// }

// static inline int compare_r(const void *a, const void *b, void *c)
// {
//     return ((int *)c)[*(int *)a] - ((int *)c)[*(int *)b];
// }

// int *run_ls_par(graph g, int *IS, int step, int nc, long long offset)
// {
//     int *C = malloc(sizeof(int) * g.N);
//     for (int i = 0; i < g.N; i++)
//         C[i] = -1;

//     long long best = 0;

// #pragma omp parallel shared(best)
//     {
//         int tid = omp_get_thread_num();
//         int nt = omp_get_num_threads();

//         local_search *ls = local_search_init(g, tid);

//         for (int i = 0; i < g.N; i++)
//             if (IS[i])
//                 local_search_add_vertex(g, ls, i);

//         local_search_greedy(g, ls);

//         int nr = 2, li = 0;
//         long long cost = ls->c;

// #pragma omp barrier

//         while (1)
//         {

// #pragma omp barrier
// #pragma omp single
//             {
//                 for (int j = 0; j < g.N; j++)
//                     C[j] = 0;
//             }
// #pragma omp barrier
// #pragma omp critical
//             {
//                 for (int j = 0; j < g.N; j++)
//                     C[j] += ls->IS[j];
//             }
// #pragma omp barrier

//             int nl = 0;
//             for (int j = 0; j < g.N; j++)
//             {
//                 if (C[j] == nt)
//                 {
//                     local_search_lock_vertex(g, ls, j);
//                     for (int _j = g.V[j]; _j < g.V[j + 1]; _j++)
//                         local_search_lock_vertex(g, ls, g.E[_j]);
//                     nl += g.V[j + 1] - g.V[j] + 1;
//                 }
//             }

//             local_search_explore(g, ls, step, g.N / nc, 0);

// #pragma omp barrier

//             for (int j = 0; j < g.N; j++)
//                 ls->T[j] = 0;

//             local_search_explore(g, ls, step, g.N / nc, 0);

//             if (ls->c > cost)
//             {
//                 cost = ls->c;
//                 li = 0;
//             }
//             else
//             {
//                 li++;
//                 if (li == 16)
//                 {
//                     li = 0;
//                     nc--;
//                     if (nc < 1)
//                         nc = 1;
//                 }
//             }

// #pragma omp barrier

//             long long validated_size = mwis_validate(g, ls->IS);

// #pragma omp critical
//             {

//                 if (validated_size > best)
//                     best = validated_size;

//                 if (tid == 0)
//                 {
//                     printf("\r%lld %d  ", best + offset, nc);
//                     fflush(stdout);
//                 }
//             }

// #pragma omp barrier
//         }

//         local_search_free(ls);
//     }
//     free(C);
//     return NULL;
// }

// int *run_ls(graph g, int *IS, int step, int nc)
// {
//     local_search *ls = local_search_init(g, 0);

//     for (int i = 0; i < g.N; i++)
//         if (IS[i])
//             local_search_add_vertex(g, ls, i);

//     local_search_greedy(g, ls);

//     int nr = 2, li = 0;
//     long long cost = ls->c;
//     while (1)
//     {
//         local_search_explore(g, ls, step, nr, 0);

//         if (ls->c > cost)
//         {
//             printf("%lld %d %d\n", ls->c, nr, li);
//             cost = ls->c;
//             li = 0;
//         }
//         else
//         {
//             li++;
//             if (li == nc)
//             {
//                 nr++;
//                 li = 0;
//             }
//         }
//     }

//     int *S = malloc(sizeof(int) * g.N);
//     for (int i = 0; i < g.N; i++)
//         S[i] = ls->IS[i];

//     local_search_free(ls);

//     return S;
// }

#include "clique_graph.h"
#include "clique_local_search.h"

// int main(int argc, char **argv)
// {
//     FILE *f = fopen(argv[1], "r");
//     clique_graph *g = clique_graph_parse_json(f);
//     fclose(f);

//     int valid = clique_graph_validate(g);

//     printf("%s Valid=%d, N=%d, M=%d, C=%d, avg D=%d\n", argv[1], valid, g->N, g->V[g->N], g->NC, g->V[g->N] / g->N);

//     clique_local_search *ls = clique_local_search_init(g, 0);

//     clique_local_search_in_order_solution(g, ls);
//     clique_local_search_greedy(g, ls);
//     clique_local_search_explore(g, ls, 100000000, 1);

//     clique_local_search_free(ls);
//     clique_graph_free(g);

//     return 0;
// }

int main(int argc, char **argv)
{
    FILE *f = fopen(argv[1], "r");
    graph *g = graph_parse(f);
    fclose(f);

    int valid = graph_validate(g);
    if (!valid)
    {
        printf("Error in graph\n");
        graph_free(g);
        return 0;
    }

    int p = 0;
    for (int i = 0; argv[1][i] != '\0'; i++)
        if (argv[1][i] == '/')
            p = i + 1;

    local_search *ls = local_search_init(g, 0);

    local_search_in_order_solution(g, ls);
    local_search_greedy(g, ls);

    local_search_explore(g, ls, 360, 0);
    long long cost_1 = mwis_validate(g, ls->independent_set);
    local_search_explore(g, ls, 1440, 0);
    long long cost_2 = mwis_validate(g, ls->independent_set);
    local_search_explore(g, ls, 1800, 0);
    long long cost_3 = mwis_validate(g, ls->independent_set);

    printf("%s %d %d %lld %lld %lld\n", argv[1] + p, g->N, g->V[g->N] / 2, cost_1, cost_2, cost_3);

    local_search_free(ls);
    graph_free(g);

    return 0;
}

// int main(int argc, char **argv)
// {
//     FILE *f = fopen(argv[1], "r");
//     graph g = graph_parse(f);
//     fclose(f);

//     int *IS = malloc(sizeof(int) * g.N);
//     int *A = malloc(sizeof(int) * g.N);
//     long long offset = 0;

//     for (int i = 0; i < g.N; i++)
//     {
//         IS[i] = 0;
//         A[i] = 1;
//     }

//     long long is = mwis_validate(g, IS);
//     int valid = 1; // graph_validate(g.N, g.V, g.E);
//     printf("%s valid graph %d |V|=%d, |E|=%d IS=%lld\n", argv[1], valid, g.N, g.V[g.N] / 2, is);

//     // kernelize_csr(g.N, g.V, g.E, g.W, A, IS, &offset, 2,
//     //               reduction_neighborhood_csr,
//     //               reduction_unconfined_csr);

//     // int *rm = malloc(sizeof(int) * g.N);
//     // graph kernel = graph_subgraph(g, A, rm);

//     // printf("Kernel %d/%d %lld (%.2f%%)\n", kernel.N, g.N, offset, ((double)(g.N - kernel.N) / (double)g.N) * 100.0);

//     // run_ls(g, IS, g.N, 10);
//     run_ls_par(g, IS, g.N / 128, 128, 0);

//     // local_search *ls = local_search_init(g, 0);
//     // local_search_explore(g, ls, 100000000, g.N / 100, 1);

//     free(IS);
//     free(A);
//     // free(rm);

//     graph_free(g);
//     // graph_free(kernel);

//     return 0;
// }

// if (argc > 4)
// {
//     f = fopen(argv[4], "r");
//     int u;
//     while (fscanf(f, "%d\n", &u) == 1)
//         IS[u - 1] = 1;
//     fclose(f);
// }