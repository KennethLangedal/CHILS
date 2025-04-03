#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <omp.h>

#include "graph.h"
#include "local_search.h"
#include "chils.h"

long long mwis_validate(graph *g, int *independent_set)
{
    long long cost = 0;
    for (int u = 0; u < g->n; u++)
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

const char *help = "CHILS --- Concurrent Hybrid Iterated Local Search\n"
                   "\nThe output of the program without -v or -b is a single line on the form:\n"
                   "instance_name,#vertices,#edges,is_weight,solution_time,total_time\n"
                   "\nThe output with -b also includes best solution after 10\%, 50\%, and 100\% of the max time/iterations\n"
                   "instance_name,#vertices,#edges,is_weight_after_10\%,is_weight_after_50\%,is_weight_after_100\%,solution_time,total_time\n"
                   "\n-h \t\tDisplay this help message\n"
                   "-v \t\tVerbose mode, output continous updates to STDOUT\n"
                   "-b \t\tBlocked mode, output additional results after 10\%, 50\%, and 100\% of the  max time/iterations\n"
                   "-g path* \tPath to the input graph in METIS format\n"
                   "-i path \tPath to initial solution (1-indexed list)\n"
                   "-o path \tPath to store the best solution found \t\t default not stored\n"
                   "-p N \t\tRun CHILS with N concurrent solutions \t\t default 1 (only local search)\n"
                   "-t sec \t\tTimout in seconds \t\t\t\t default 3600 seconds\n"
                   "-s sec \t\tAlternating interval for CHILS \t\t\t default 5 seconds\n"
                   "-q N \t\tMax queue size after perturbe \t\t\t default 32\n"
                   "-c T \t\tSet a specific number of threads  \t\t default OMP_NUM_THREADS\n"
                   "-r s \t\tSet a specific random seed \t\t\t default time(NULL)\n"
                   "\n"
                   "-n it \t\tMax CHILS iterations \t\t\t\t default inf\n"
                   "-m it \t\tMax local search iterations \t\t\t default inf\n"
                   "\n* Mandatory input";

int main(int argc, char **argv)
{
    char *graph_path = NULL,
         *initial_solution_path = NULL,
         *solution_path = NULL;
    int verbose = 0, blocked = 0, run_chils = 1, max_queue = 32, num_threads = -1;
    double timeout = 3600, step = 5, reduction_timout = 30;

    long long cl = LLONG_MAX, il = LLONG_MAX;

    unsigned int seed = time(NULL);

    int command;

    while ((command = getopt(argc, argv, "hvbg:i:o:p:t:n:s:m:q:c:r:")) != -1)
    {
        switch (command)
        {
        case 'h':
            printf("%s\n", help);
            return 0;
        case 'v':
            verbose = 1;
            break;
        case 'b':
            blocked = 1;
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
            run_chils = atoi(optarg);
            break;
        case 't':
            timeout = atof(optarg);
            break;
        case 'n':
            cl = atoll(optarg);
            break;
        case 's':
            step = atof(optarg);
            break;
        case 'm':
            il = atoll(optarg);
            break;
        case 'q':
            max_queue = atoi(optarg);
            break;
        case 'c':
            num_threads = atoi(optarg);
            break;
        case 'r':
            seed = atoi(optarg);
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

        initial_solution = malloc(sizeof(int) * g->n);
        for (int i = 0; i < g->n; i++)
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
        if (run_chils <= 1)
            printf("Running iterated local search\n");
        else
            printf("Running CHILS using %d concurrent solutions\n", run_chils);

        if (blocked)
            printf("Blocked mode\n");

        printf("Input: \t\t\t%s\n", graph_path + path_offset);
        printf("Vertices: \t\t%d\n", g->n);
        printf("Edges: \t\t\t%d\n", g->V[g->n] / 2);
        printf("Seed: \t\t\t%u\n", seed);
        if (solution_path != NULL)
            printf("Output: \t\t%s\n", solution_path);
        if (cl < LLONG_MAX)
            printf("Tiomeout: \t\t%.2lf seconds or %lld CHILS iterations\n", timeout, cl);
        else
            printf("Tiomeout: \t\t%.2lf seconds\n", timeout);
        printf("Max queue size: \t%d\n", max_queue);
        if (run_chils > 1)
        {
            if (il < LLONG_MAX)
                printf("CHILS interval: \t%.2lf seconds or %lld iterations\n", step, il);
            else
                printf("CHILS interval: \t%.2lf seconds\n", step);
        }
        else if (il < LLONG_MAX)
        {
            printf("LS iteration limit: \t%lld\n", il);
        }
        if (initial_solution_path != NULL)
            printf("Initial solution: \t%lld\n", initial_solution_weight);
        printf("\n");
    }

    long long w10, w50, w100;
    double t10 = timeout * 0.1, t50 = timeout * 0.4, t100 = timeout * 0.5, tb = 0.0, t_total = 0.0;
    long long c10 = cl / 10ll, c50 = (cl / 10ll) * 4, c100 = (cl / 10ll) * 5;
    long long i10 = il / 10ll, i50 = (il / 10ll) * 4, i100 = (il / 10ll) * 5;

    int *solution = malloc(sizeof(int) * g->n);
    for (int i = 0; i < g->n; i++)
        solution[i] = 0;

    if (run_chils > 1)
    {
        if (num_threads > 0)
            omp_set_num_threads(num_threads);

        int nt = omp_get_max_threads();
        printf("%d\n", nt);
        if (nt > run_chils)
            omp_set_num_threads(run_chils);

        chils *c = chils_init(g, run_chils, seed);
        c->step_time = step;
        c->step_count = il;

        if (initial_solution != NULL)
            chils_set_solution(g, c, initial_solution);

        for (int i = 0; i < run_chils; i++)
        {
            c->LS[i]->max_queue = max_queue + (4 * i);
            c->LS_core[i]->max_queue = max_queue + (4 * i);
        }

        double start = omp_get_wtime();

        if (blocked)
        {
            chils_run(g, c, t10, c10, verbose);
            w10 = mwis_validate(g, chils_get_best_independent_set(c));
            chils_run(g, c, t50, c50, verbose);
            w50 = mwis_validate(g, chils_get_best_independent_set(c));
            chils_run(g, c, t100, c100, verbose);
            w100 = mwis_validate(g, chils_get_best_independent_set(c));
        }
        else
        {
            chils_run(g, c, timeout, cl, verbose);
            w100 = mwis_validate(g, chils_get_best_independent_set(c));
        }

        double end = omp_get_wtime();
        t_total = end - start;

        tb = c->time;

        int *best = chils_get_best_independent_set(c);
        for (int i = 0; i < g->n; i++)
            solution[i] = best[i];

        chils_free(c);
    }
    else
    {
        local_search *ls = local_search_init(g, seed);

        if (initial_solution != NULL)
            for (int u = 0; u < g->n; u++)
                if (initial_solution[u])
                    local_search_add_vertex(g, ls, u);

        ls->max_queue = max_queue;

        double start = omp_get_wtime();

        if (blocked)
        {
            local_search_explore(g, ls, t10, i10, verbose);
            w10 = mwis_validate(g, ls->independent_set);
            local_search_explore(g, ls, t50, i50, verbose);
            w50 = mwis_validate(g, ls->independent_set);
            local_search_explore(g, ls, t100, i100, verbose);
            w100 = mwis_validate(g, ls->independent_set);
        }
        else
        {
            local_search_explore(g, ls, timeout, il, verbose);
            w100 = mwis_validate(g, ls->independent_set);
        }

        double end = omp_get_wtime();
        t_total = end - start;

        tb = ls->time;

        for (int i = 0; i < g->n; i++)
            solution[i] = ls->independent_set[i];

        local_search_free(ls);
    }

    if (blocked)
        printf("%s,%d,%d,%lld,%lld,%lld,%.4lf,%.4lf\n", graph_path + path_offset,
               g->n, g->V[g->n] / 2, w10, w50, w100, tb, t_total);
    else
        printf("%s,%d,%d,%lld,%.4lf,%.4lf\n", graph_path + path_offset,
               g->n, g->V[g->n] / 2, w100, tb, t_total);

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

            for (int i = 0; i < g->n; i++)
                if (solution[i])
                    fprintf(f, "%d\n", i + 1);

            fclose(f);
        }
    }

    free(solution);
    free(initial_solution);
    graph_free(g);

    return 0;
}
