#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

#include "graph.h"
#include "reductions.h"
#include "local_search.h"
#include "pils.h"

long long mwis_validate(graph *g, int *independent_set)
{
    long long cost = 0;
    for (int u = 0; u < g->N; u++)
    {
        if (!independent_set[u])
            continue;

        cost += g->W[u];
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (independent_set[v])
                return -1;
        }
    }
    return cost;
}

const char *help = "PILS --- Parallel Iterated Local Search\n"
                   "\nThe output of the program without -v is a single line on the form:\n"
                   "instance_name #vertices #edges W_after_10% W_after_50% W_after_t\n"
                   "\n-h \t\tDisplay this help message\n"
                   "-v \t\tVerbose mode, output continous updates to STDOUT\n"
                   "-g path* \tPath to the input graph in METIS format\n"
                   "-i path* \tPath to initial solution (1-indexed list)\n"
                   "-o path \tPath to store the best solution found \t\t default not stored\n"
                   "-p N \t\tRun PILS with N concurrent solutions \t\t default 1 (only local search)\n"
                   "-t sec \t\tTimout in seconds \t\t\t\t default 3600 seconds\n"
                   "-s sec \t\tAlternating interval for PILS \t\t\t default 5 seconds\n"
                   "-r sec \t\tTime to spend on initial reductions \t\t default 0 seconds\n"
                   "-q N \t\tMax queue size after shake \t\t\t default 32\n"
                   "-c T \t\tSet a specific number of threads  \t\t default OMP_NUM_THREADS\n"
                   "\n* Mandatory input";

int main(int argc, char **argv)
{
    char *graph_path = NULL,
         *initial_solution_path = NULL,
         *solution_path = NULL;
    int verbose = 0, run_pils = 1, run_reductions = 0, max_queue = 32, num_threads = -1;
    double timeout = 3600, step = 5, reduction_timout = 30;

    int command;

    while ((command = getopt(argc, argv, "hvg:i:o:p:t:s:r:q:c:")) != -1)
    {
        switch (command)
        {
        case 'h':
            printf("%s\n", help);
            return 0;
        case 'v':
            verbose = 1;
            break;
        case 'g':
            graph_path = optarg;
            break;
        case 'i':
            initial_solution_path = optarg;
            break;
        case 'o':
            solution_path = optarg;
            break;
        case 'p':
            run_pils = atoi(optarg);
            break;
        case 't':
            timeout = atof(optarg);
            break;
        case 's':
            step = atof(optarg);
            break;
        case 'r':
            run_reductions = 1;
            reduction_timout = atof(optarg);
            break;
        case 'q':
            max_queue = atoi(optarg);
            break;
        case 'c':
            num_threads = atoi(optarg);
            break;
        case '?':
            return 1;

        default:
            return 1;
        }
    }

    if (graph_path == NULL)
    {
        fprintf(stderr, "Input graph is required, run with -h for more information\n");
        return 1;
    }

    FILE *f = fopen(graph_path, "r");
    if (f == NULL)
    {
        fprintf(stderr, "Unable to open file %s\n", graph_path);
        return 1;
    }
    graph *g = graph_parse(f);
    fclose(f);

    if (!graph_validate(g))
    {
        fprintf(stderr, "Errors in input graph\n");
        return 1;
    }

    int *initial_solution = NULL;
    long long initial_solution_weight = 0;

    if (initial_solution_path != NULL)
    {
        f = fopen(initial_solution_path, "r");
        if (f == NULL)
        {
            fprintf(stderr, "Unable to open file %s\n", initial_solution_path);
            return 1;
        }

        initial_solution = malloc(sizeof(int) * g->N);
        for (int i = 0; i < g->N; i++)
            initial_solution[i] = 0;

        int u = 0;
        while (fscanf(f, "%d", &u) != EOF)
            initial_solution[u - 1] = 1;

        initial_solution_weight = mwis_validate(g, initial_solution);
        if (initial_solution_weight < 0)
        {
            fprintf(stderr, "Initial solution is not valid\n");
            return 1;
        }

        fclose(f);
    }

    int path_offset = 0, path_end = 0;
    for (int i = 0; graph_path[i] != '\0'; i++)
        if (graph_path[i] == '/')
            path_offset = i + 1;
        else if (graph_path[i] == '.')
            path_end = i;
    graph_path[path_end] = '\0';

    if (verbose)
    {
        if (run_pils <= 1)
            printf("Running iterated local search\n");
        else
            printf("Running PILS using %d concurrent solutions\n", run_pils);

        printf("Input: \t\t\t%s\n", graph_path + path_offset);
        printf("Vertices: \t\t%d\n", g->N);
        printf("Edges: \t\t\t%d\n", g->V[g->N] / 2);
        if (solution_path != NULL)
            printf("Output: \t\t%s\n", solution_path);
        printf("Tiomeout: \t\t%.2lf seconds\n", timeout);
        printf("Max queue size: \t%d\n", max_queue);
        if (run_pils > 1)
            printf("PILS interval: \t\t%.2lf seconds\n", step);
        if (initial_solution_path != NULL)
            printf("Initial solution: \t%lld\n", initial_solution_weight);
    }

    int *reverse_mapping = NULL;
    int *original_solution = NULL;
    graph *original_graph = NULL;
    long long offset = 0;

    if (run_reductions && initial_solution == NULL)
    {
        reverse_mapping = malloc(sizeof(int) * g->N);
        original_solution = malloc(sizeof(int) * g->N);
        for (int i = 0; i < g->N; i++)
            original_solution[i] = 0;
        int *A = malloc(sizeof(int) * g->N);
        for (int i = 0; i < g->N; i++)
            A[i] = 1;

        kernelize_csr(g->N, g->V, g->E, g->W, A,
                      original_solution, &offset, reduction_timout, 7,
                      reduction_neighborhood_csr,
                      reduction_clique_csr,
                      reduction_domination_csr,
                      reduction_single_edge_csr,
                      reduction_extended_single_edge_csr,
                      reduction_extended_twin_csr,
                      reduction_extended_unconfined_csr);

        original_graph = g;
        g = graph_subgraph(original_graph, A, reverse_mapping);

        if (verbose)
            printf("Kernel: \t\t|V|=%d |E|=%d\n", g->N, g->V[g->N] / 2);

        free(A);
    }

    long long w10, w50, w100;
    double t10 = timeout * 0.1, t50 = timeout * 0.4, t100 = timeout * 0.5;

    int *solution = malloc(sizeof(int) * g->N);

    if (run_pils > 1)
    {
        if (num_threads > 0)
            omp_set_num_threads(num_threads);

        pils *p = pils_init(g, run_pils);
        p->step = step;

        if (initial_solution != NULL)
            pils_set_solution(g, p, initial_solution);

        for (int i = 0; i < run_pils; i++)
            p->LS[i]->max_queue = max_queue + (i * 8);

        pils_run(g, p, t10, verbose, offset);
        w10 = mwis_validate(g, pils_get_best_independent_set(p)) + offset;
        pils_run(g, p, t50, verbose, offset);
        w50 = mwis_validate(g, pils_get_best_independent_set(p)) + offset;
        pils_run(g, p, t100, verbose, offset);
        w100 = mwis_validate(g, pils_get_best_independent_set(p)) + offset;

        int *best = pils_get_best_independent_set(p);
        for (int i = 0; i < g->N; i++)
            solution[i] = best[i];

        pils_free(p);
    }
    else
    {
        local_search *ls = local_search_init(g, 0);

        if (initial_solution != NULL)
            for (int u = 0; u < g->N; u++)
                if (initial_solution[u])
                    local_search_add_vertex(g, ls, u);

        ls->max_queue = max_queue;
        local_search_explore(g, ls, t10, verbose, offset);
        w10 = mwis_validate(g, ls->independent_set) + offset;
        local_search_explore(g, ls, t50, verbose, offset);
        w50 = mwis_validate(g, ls->independent_set) + offset;
        local_search_explore(g, ls, t100, verbose, offset);
        w100 = mwis_validate(g, ls->independent_set) + offset;

        for (int i = 0; i < g->N; i++)
            solution[i] = ls->independent_set[i];

        local_search_free(ls);
    }

    printf("%s,%d,%d,%lld,%lld,%lld\n", graph_path + path_offset,
           g->N, g->V[g->N] / 2, w10, w50, w100);

    if (solution_path != NULL)
    {
        if (original_solution != NULL)
        {
            for (int u = 0; u < g->N; u++)
                if (solution[u])
                    original_solution[reverse_mapping[u]] = 1;
        }

        f = fopen(solution_path, "w");
        if (f == NULL)
        {
            fprintf(stderr, "Unable to open file %s\n", solution_path);
        }
        else
        {
            if (verbose)
                printf("\nStoring solution of size %lld to %s\n", w100, solution_path);

            if (original_solution == NULL)
            {
                for (int i = 0; i < g->N; i++)
                    if (solution[i])
                        fprintf(f, "%d\n", i + 1);
            }
            else
            {
                for (int i = 0; i < original_graph->N; i++)
                    if (original_solution[i])
                        fprintf(f, "%d\n", i + 1);
            }

            fclose(f);
        }
    }

    free(solution);
    free(initial_solution);
    free(original_solution);
    free(reverse_mapping);
    graph_free(g);
    graph_free(original_graph);

    return 0;
}
