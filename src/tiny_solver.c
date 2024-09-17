#include "tiny_solver.h"

#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <limits.h>

#define MAX_NODES 128

tiny_solver *tiny_solver_init(int N)
{
    tiny_solver *solver = malloc(sizeof(tiny_solver));

    solver->subgraph = malloc(sizeof(char *) * MAX_NODES);
    char *data = aligned_alloc(32, sizeof(char) * MAX_NODES * MAX_NODES);
    for (int i = 0; i < MAX_NODES; i++)
        solver->subgraph[i] = data + i * MAX_NODES;

    solver->subgraph_W = aligned_alloc(32, sizeof(long long) * MAX_NODES);

    solver->forward_map = malloc(sizeof(int) * N);
    for (int i = 0; i < N; i++)
        solver->forward_map[i] = 0;
    solver->reverse_map = malloc(sizeof(int) * MAX_NODES);
    for (int i = 0; i < MAX_NODES; i++)
        solver->reverse_map[i] = 0;

    solver->independent_set = malloc(sizeof(int) * MAX_NODES);

    char *solution_data = aligned_alloc(32, sizeof(char) * MAX_NODES * MAX_NODES);
    solver->per_layer_solution = malloc(sizeof(char *) * MAX_NODES);
    for (int i = 0; i < MAX_NODES; i++)
        solver->per_layer_solution[i] = solution_data + i * MAX_NODES;
    solver->per_layer_weight = aligned_alloc(32, sizeof(long long) * MAX_NODES);
    solver->branch = malloc(sizeof(int) * MAX_NODES);

    return solver;
}

void tiny_solver_free(tiny_solver *solver)
{
    free(*solver->subgraph);
    free(solver->subgraph);
    free(solver->subgraph_W);
    free(solver->forward_map);
    free(solver->reverse_map);
    free(solver->independent_set);
    free(*solver->per_layer_solution);
    free(solver->per_layer_solution);
    free(solver->per_layer_weight);
    free(solver->branch);
    free(solver);
}

void tiny_solver_clear(tiny_solver *solver)
{
    solver->subgraph_N = 0;
    memset(*solver->subgraph, 0, sizeof(char) * MAX_NODES * MAX_NODES);

    solver->independent_set_weight = 0;
    for (int i = 0; i < MAX_NODES; i++)
        solver->independent_set[i] = 0;

    solver->time_limit_exceeded = 0;
    solver->weight_limit_exceeded = 0;
    solver->node_limit_exceeded = 0;
}

static inline void tiny_solver_reduce(tiny_solver *solver, int layer)
{
    char *S = solver->per_layer_solution[layer];
    char **g = solver->subgraph;
    long long *W = solver->subgraph_W;
    int N = solver->subgraph_N;
    int imp = 1;

    while (imp)
    {
        imp = 0;

        // Neighbourhood rule
        for (int u = 0; u < N; u++)
        {
            if (S[u] != 0)
                continue;

            long long nw = 0;
            for (int v = 0; v < N; v++)
                if (v != u && g[u][v] && S[v] == 0)
                    nw += W[v];

            if (nw <= W[u])
            {
                imp = 1;
                solver->per_layer_weight[layer] += W[u];
                S[u] = 1;
                for (int v = 0; v < N; v++)
                    if (u != v && g[u][v])
                        S[v] = -1;
            }
        }
        if (imp)
            continue;

        // Domination rule
        for (int u = 0; u < N; u++)
        {
            if (S[u] != 0)
                continue;
            for (int v = 0; v < N; v++)
            {
                if (S[v] != 0 || v == u || W[u] > W[v])
                    continue;

                int dom = 1;
                for (int i = 0; i < N; i++)
                {
                    if (i == u || i == v)
                        continue;
                    if (g[v][i] && !g[u][i])
                    {
                        dom = 0;
                        break;
                    }
                }

                if (dom)
                {
                    imp = 1;
                    S[u] = -1;
                    break;
                }
            }
        }
    }
}

void tiny_solver_solve(tiny_solver *solver, double tl, long long Wl)
{
    double start = omp_get_wtime();
    int layer = 0, vertex = 0, branch_count = 0;

    for (int i = 0; i < solver->subgraph_N; i++)
        solver->per_layer_solution[layer][i] = 0;

    solver->per_layer_weight[layer] = 0;
    solver->independent_set_weight = 0;

    solver->branch[layer] = 1;

    tiny_solver_reduce(solver, layer);

    while (layer >= 0)
    {
        if ((branch_count++ & 1023) == 0)
        {
            branch_count = 0;
            double elapsed = omp_get_wtime() - start;
            if (elapsed > tl)
            {
                solver->time_limit_exceeded = 1;
                return;
            }
        }

        // Find max degree vertex
        vertex = solver->subgraph_N;
        int max_degree = 0;
        for (int i = 0; i < solver->subgraph_N; i++)
        {
            if (solver->per_layer_solution[layer][i] != 0)
                continue;
            int degree = 0;
            for (int j = 0; j < solver->subgraph_N; j++)
                if (i != j && solver->subgraph[i][j] &&
                    solver->per_layer_solution[layer][j] == 0)
                    degree++;
            if (degree > max_degree)
            {
                max_degree = degree;
                vertex = i;
            }
        }

        // Fast upperbound
        long long ub = solver->per_layer_weight[layer];
        for (int i = 0; i < solver->subgraph_N; i++)
            if (solver->per_layer_solution[layer][i] == 0)
                ub += solver->subgraph_W[i];

        // Base case
        if (vertex == solver->subgraph_N || ub <= solver->independent_set_weight)
        {
            // New best solution
            if (solver->per_layer_weight[layer] > solver->independent_set_weight)
            {
                solver->independent_set_weight = solver->per_layer_weight[layer];
                for (int i = 0; i < solver->subgraph_N; i++)
                    solver->independent_set[i] = solver->per_layer_solution[layer][i];

                if (solver->independent_set_weight > Wl)
                {
                    solver->weight_limit_exceeded = 1;
                    return;
                }
            }

            // Backtrack to last branch
            layer--;
            while (layer >= 0 && solver->branch[layer] == 0)
                layer--;

            if (layer >= 0)
                solver->branch[layer] = 0;
        }
        // Include branch
        else if (solver->branch[layer] == 1)
        {
            layer++;
            solver->branch[layer] = 1;
            for (int i = 0; i < solver->subgraph_N; i++)
                solver->per_layer_solution[layer][i] = solver->per_layer_solution[layer - 1][i];

            solver->per_layer_solution[layer][vertex] = 1;
            solver->per_layer_weight[layer] = solver->per_layer_weight[layer - 1] + solver->subgraph_W[vertex];

            for (int u = 0; u < solver->subgraph_N; u++)
                if (u != vertex && solver->subgraph[vertex][u])
                    solver->per_layer_solution[layer][u] = -1;

            tiny_solver_reduce(solver, layer);
        }
        // Exclude branch
        else if (solver->branch[layer] == 0)
        {
            layer++;
            solver->branch[layer] = 1;
            for (int i = 0; i < solver->subgraph_N; i++)
                solver->per_layer_solution[layer][i] = solver->per_layer_solution[layer - 1][i];

            solver->per_layer_solution[layer][vertex] = -1;
            solver->per_layer_weight[layer] = solver->per_layer_weight[layer - 1];

            tiny_solver_reduce(solver, layer);
        }
    }
}

void tiny_solver_solve_neighbourhood(tiny_solver *solver, double tl, long long Wl, int u,
                                     const int *V, const int *E, const long long *W, const int *A)
{
    tiny_solver_clear(solver);

    int N = 0;
    for (int i = V[u]; i < V[u + 1]; i++)
    {
        int v = E[i];
        if (!A[v])
            continue;
        solver->forward_map[v] = N;
        solver->reverse_map[N] = v;
        solver->subgraph_W[N] = W[v];
        N++;

        if (N > MAX_NODES)
        {
            solver->node_limit_exceeded = 1;
            return;
        }
    }

    solver->subgraph_N = N;

    for (int i = V[u]; i < V[u + 1]; i++)
    {
        int v = E[i];
        if (!A[v])
            continue;

        solver->subgraph[solver->forward_map[u]][solver->forward_map[u]] = 1;

        for (int j = V[v]; j < V[v + 1]; j++)
        {
            int w = E[j];
            if (!A[w] || solver->forward_map[w] < 0 || solver->forward_map[w] >= N ||
                solver->reverse_map[solver->forward_map[w]] != w)
                continue;

            solver->subgraph[solver->forward_map[v]][solver->forward_map[w]] = 1;
        }
    }

    tiny_solver_solve(solver, tl, Wl);
}

void tiny_solver_build_subgraph(tiny_solver *solver, int sN, const int *sV,
                                const int *V, const int *E, const long long *W, const int *A)
{
    tiny_solver_clear(solver);

    int N = 0;
    for (int i = 0; i < sN; i++)
    {
        int u = sV[i];
        if (!A[u])
            continue;
        solver->forward_map[u] = N;
        solver->reverse_map[N] = u;
        solver->subgraph_W[N] = W[u];
        N++;

        if (N > MAX_NODES)
        {
            solver->node_limit_exceeded = 1;
            return;
        }
    }

    solver->subgraph_N = N;

    for (int i = 0; i < sN; i++)
    {
        int u = sV[i];
        if (!A[u])
            continue;

        solver->subgraph[solver->forward_map[u]][solver->forward_map[u]] = 1;

        for (int j = V[u]; j < V[u + 1]; j++)
        {
            int v = E[j];
            if (!A[v] || solver->forward_map[v] < 0 || solver->forward_map[v] >= N ||
                solver->reverse_map[solver->forward_map[v]] != v)
                continue;

            solver->subgraph[solver->forward_map[u]][solver->forward_map[v]] = 1;
        }
    }
}

void tiny_solver_solve_subgraph(tiny_solver *solver, double tl, long long Wl, int sN, const int *sV,
                                const int *V, const int *E, const long long *W, const int *A)
{
    tiny_solver_build_subgraph(solver, sN, sV, V, E, W, A);
    tiny_solver_solve(solver, tl, Wl);
}