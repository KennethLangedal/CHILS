#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

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
                   "-o path \tPath to store the best solution found \t\t default not stored\n"
                   "-p N \t\tRun PILS with N concurrent solutions \t\t default 1 (only local search)\n"
                   "-t sec \t\tTimout in seconds \t\t\t\t default 3600 seconds\n"
                   "-s sec \t\tAlternating interval for PILS \t\t\t default 10 seconds\n"
                   "\n* Mandatory input";

int main(int argc, char **argv)
{
    char *graph_path = NULL, *solution_path = NULL;
    int verbose = 0, run_pils = 0;
    double timeout = 3600, step = 10;

    int command;

    while ((command = getopt(argc, argv, "hvg:o:p:t:s:")) != -1)
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
        graph_free(g);
        return 1;
    }

    int path_offset = 0;
    for (int i = 0; graph_path[i] != '\0'; i++)
        if (graph_path[i] == '/')
            path_offset = i + 1;

    if (verbose)
    {
        if (run_pils == 0)
            printf("Running iterated local search\n");
        else
            printf("Running PILS using %d concurrent solutions\n", run_pils);

        printf("Input: \t\t%s\n", graph_path + path_offset);
        printf("Vertices: \t%d\n", g->N);
        printf("Edges: \t\t%d\n", g->V[g->N] / 2);
        if (solution_path != NULL)
            printf("Output: \t%s\n", solution_path);
        printf("Tiomeout: \t%.2lf seconds\n", timeout);
        if (run_pils > 0)
            printf("PILS interval: \t%.2lf seconds\n", step);
        printf("\n");
    }

    long long w10, w50, w100;
    double t10 = timeout * 0.1, t50 = timeout * 0.4, t100 = timeout * 0.5;

    int *solution = malloc(sizeof(int) * g->N);

    if (run_pils > 0)
    {
        pils *p = pils_init(g, run_pils);
        p->step_full = step;
        p->step_reduced = step;

        pils_run(g, p, t10, verbose);
        w10 = mwis_validate(g, pils_get_best_independent_set(p));
        pils_run(g, p, t50, verbose);
        w50 = mwis_validate(g, pils_get_best_independent_set(p));
        pils_run(g, p, t100, verbose);
        w100 = mwis_validate(g, pils_get_best_independent_set(p));

        int *best = pils_get_best_independent_set(p);
        for (int i = 0; i < g->N; i++)
            solution[i] = best[i];

        pils_free(p);
    }
    else
    {
        local_search *ls = local_search_init(g, time(NULL));

        local_search_explore(g, ls, t10, verbose);
        w10 = mwis_validate(g, ls->independent_set);
        local_search_explore(g, ls, t50, verbose);
        w50 = mwis_validate(g, ls->independent_set);
        local_search_explore(g, ls, t100, verbose);
        w100 = mwis_validate(g, ls->independent_set);

        for (int i = 0; i < g->N; i++)
            solution[i] = ls->independent_set[i];

        local_search_free(ls);
    }

    if (!verbose)
        printf("%s %d %d %lld %lld %lld\n", graph_path + path_offset,
               g->N, g->V[g->N] / 2, w10, w50, w100);

    if (solution_path != NULL)
    {
        f = fopen(solution_path, "w");
        if (f == NULL)
        {
            fprintf(stderr, "Unable to open file %s\n", solution_path);
        }
        else
        {
            if (verbose)
                printf("\nStoring solution of size %lld to %s\n", w100, solution_path);

            for (int i = 0; i < g->N; i++)
                if (solution[i])
                    fprintf(f, "%d\n", i + 1);
            fclose(f);
        }
    }

    free(solution);
    graph_free(g);

    return 0;
}
