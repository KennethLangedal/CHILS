#pragma once

typedef struct
{
    // Subgraph
    int subgraph_N;
    char **subgraph;
    long long *subgraph_W;

    // Vertex mapping
    int *forward_map, *reverse_map;

    // Solution
    int *independent_set;
    long long independent_set_weight;

    // Flags
    int time_limit_exceeded,
        weight_limit_exceeded,
        node_limit_exceeded;

    // Internal structures
    char **per_layer_solution;
    long long *per_layer_weight;
    int *branch;
} tiny_solver;

tiny_solver *tiny_solver_init(int N);

void tiny_solver_free(tiny_solver *solver);

void tiny_solver_clear(tiny_solver *solver);

void tiny_solver_solve(tiny_solver *solver, double tl, long long Wl);

void tiny_solver_solve_neighbourhood(tiny_solver *solver, double tl, long long Wl, int u,
                                     const int *V, const int *E, const long long *W, const int *A);

void tiny_solver_build_subgraph(tiny_solver *solver, int sN, const int *sV,
                                const int *V, const int *E, const long long *W, const int *A);

void tiny_solver_solve_subgraph(tiny_solver *solver, double tl, long long Wl, int sN, const int *sV,
                                const int *V, const int *E, const long long *W, const int *A);
