#include "graph.h"

#include <omp.h>
#include <stdlib.h>
#include <limits.h>

static inline void parse_id(char *Data, size_t *p, long long *v)
{
    while (Data[*p] < '0' || Data[*p] > '9')
        (*p)++;

    *v = 0;
    while (Data[*p] >= '0' && Data[*p] <= '9')
        *v = (*v) * 10 + Data[(*p)++] - '0';
}

static inline void skip_line(char *Data, size_t *p)
{
    while (Data[*p] != '\n')
        (*p)++;
    (*p)++;
}

static inline int graph_compare(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

graph *graph_parse(FILE *f)
{
    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *Data = malloc(size);
    size_t red = fread(Data, 1, size, f);
    size_t p = 0;

    while (Data[p] == '%')
        skip_line(Data, &p);

    long long n, m, t = 0;
    parse_id(Data, &p, &n);
    parse_id(Data, &p, &m);
    while (Data[p] == ' ')
        p++;
    if (Data[p] >= '0' && Data[p] <= '9')
        parse_id(Data, &p, &t);

    skip_line(Data, &p);

    int vertex_weights = t >= 10,
        edge_weights = (t == 1 || t == 11);

    if (n >= INT_MAX)
    {
        fprintf(stderr, "Number of vertices must be less than %d, got %lld\n", INT_MAX, n);
        exit(1);
    }

    long long *V = malloc(sizeof(long long) * (n + 1));
    int *E = malloc(sizeof(int) * (m * 2));

    long long *W = malloc(sizeof(long long) * n);

    long long ei = 0;
    for (int u = 0; u < n; u++)
    {
        W[u] = 1;
        V[u] = ei;

        while (p < size && Data[p] == '%')
            skip_line(Data, &p);

        if (vertex_weights)
            parse_id(Data, &p, W + u);

        while (ei < m * 2)
        {
            while (Data[p] == ' ')
                p++;
            if (Data[p] == '\n' || Data[p] == EOF)
                break;

            long long e;
            parse_id(Data, &p, &e);

            if (e > n || e <= 0)
            {
                fprintf(stderr, "Edge endpoint out of bounds, {%lld, %lld}\n", u + 1ll, e);
                exit(1);
            }

            E[ei++] = e - 1;

            if (edge_weights)
                parse_id(Data, &p, &e);
        }
        p++;

        qsort(E + V[u], ei - V[u], sizeof(int), graph_compare);
    }
    V[n] = ei;

    free(Data);

    graph *g = malloc(sizeof(graph));
    *g = (graph){.n = n, .m = m * 2, .V = V, .E = E, .W = W};

    return g;
}

// store graph in metis format
void graph_store(FILE *f, graph *g)
{
    fprintf(f, "%d %lld 10\n", g->n, g->m / 2);
    for (int u = 0; u < g->n; u++)
    {
        fprintf(f, "%lld", g->W[u]);
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
            fprintf(f, " %d", g->E[i] + 1);
        fprintf(f, "\n");
    }
}

void graph_free(graph *g)
{
    if (g == NULL)
        return;

    free(g->V);
    free(g->E);
    free(g->W);

    free(g);
}

int graph_validate(graph *g)
{
    long long *Edge_pos = malloc(sizeof(long long) * g->n);

    long long M = 0;
    for (int u = 0; u < g->n; u++)
    {
        Edge_pos[u] = g->V[u + 1];

        long long d_u = g->V[u + 1] - g->V[u];
        if (d_u < 0)
        {
            fprintf(stderr, "Error in neighborhood list V: Vertex %d starts at position "
                            "%lld and ends at position %lld\n",
                    u + 1, g->V[u], g->V[u + 1]);
            return 0;
        }

        M += d_u;

        int first = 1;
        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
        {
            if (i < 0 || i >= g->m)
            {
                fprintf(stderr, "Error in neighborhood list V: Vertex %d starts at position "
                                "%lld and ends at position %lld\n",
                        u + 1, g->V[u], g->V[u + 1]);
                return 0;
            }

            int v = g->E[i];
            if (v < 0 || v >= g->n)
            {
                fprintf(stderr, "Edge endpoint out of bounds for {%d,%d}\n", u + 1, v + 1);
                return 0;
            }
            if (v == u)
            {
                fprintf(stderr, "Self edges are not allowd {%d,%d}\n", u + 1, u + 1);
                return 0;
            }
            if (i > g->V[u] && v <= g->E[i - 1])
            {
                fprintf(stderr, "Unsorted neighborhood for vertex %d: {...,%d,%d,...}\n", u + 1, g->E[i - 1] + 1, v + 1);
                return 0;
            }

            if (u > v)
            {
                if (Edge_pos[v] >= g->V[v + 1] || g->E[Edge_pos[v]] != u)
                {
                    if (Edge_pos[v] >= g->V[v + 1] || g->E[Edge_pos[v]] > u)
                        fprintf(stderr, "Undirected edge encountered: Found {%d,%d} but not {%d,%d}\n", u + 1, v + 1, v + 1, u + 1);
                    else
                        fprintf(stderr, "Undirected edge encountered: Found {%d,%d} but not {%d,%d}\n", v + 1, g->E[Edge_pos[v]] + 1, g->E[Edge_pos[v]] + 1, v + 1);
                    return 0;
                }
                Edge_pos[v]++;
            }
            else if (first)
            {
                Edge_pos[u] = i;
                first = 0;
            }
        }
    }

    for (int u = 0; u < g->n; u++)
    {
        if (Edge_pos[u] != g->V[u + 1])
        {
            int v = g->E[Edge_pos[u]];
            fprintf(stderr, "Undirected edge encountered: Found {%d,%d} but not {%d,%d}\n", u + 1, v + 1, v + 1, u + 1);
            return 0;
        }
    }

    if (M != g->V[g->n] || M != g->m)
    {
        fprintf(stderr, "Wrong edge count, found %lld, but file says %lld\n", M / 2, g->m / 2);
        return 0;
    }

    free(Edge_pos);

    return 1;
}

graph *graph_subgraph(graph *g, int *Mask, int *RM)
{
    int *FM = malloc(sizeof(int) * g->n);
    long long n = 0, m = 0;
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        FM[u] = n;
        RM[n] = u;
        n++;

        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
            if (Mask[g->E[i]])
                m++;
    }

    graph *sg = malloc(sizeof(graph));
    *sg = (graph){.n = n};

    sg->V = malloc(sizeof(long long) * (n + 1));
    sg->E = malloc(sizeof(int) * m);
    sg->W = malloc(sizeof(long long) * n);

    m = 0;
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        sg->W[FM[u]] = g->W[u];
        sg->V[FM[u]] = m;

        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (!Mask[v])
                continue;

            sg->E[m] = FM[v];
            m++;
        }
    }
    sg->V[sg->n] = m;

    free(FM);

    return sg;
}

void graph_subgraph_par(graph *g, graph *sg, int *Mask, int *RM, int *FM, long long *S1, long long *S2)
{
    int nt = omp_get_num_threads();
    int tid = omp_get_thread_num();
    long long n = 0, m = 0;

#pragma omp for nowait
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        n++;

        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
            if (Mask[g->E[i]])
                m++;
    }

    S1[tid] = n;
    S2[tid] = m;

#pragma omp barrier

    long long n_offset = 0;
    for (int i = 0; i < tid; i++)
        n_offset += S1[i];

    n = n_offset;
#pragma omp for
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        FM[u] = n;
        RM[n] = u;
        sg->W[n] = g->W[u];
        n++;
    }

    long long m_offset = 0;
    for (int i = 0; i < tid; i++)
        m_offset += S2[i];

    m = m_offset;
#pragma omp for nowait
    for (int u = 0; u < g->n; u++)
    {
        if (!Mask[u])
            continue;

        sg->V[FM[u]] = m;

        for (long long i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            if (!Mask[v])
                continue;

            sg->E[m] = FM[v];
            m++;
        }
    }

    if (tid == nt - 1)
    {
        sg->n = n;
        sg->m = m;
        sg->V[sg->n] = m;
    }
#pragma omp barrier
}