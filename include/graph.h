#pragma once
#include <stdio.h>

typedef struct
{
    int N;
    int *V, *E;
    long long *W;
} graph;

graph *graph_parse(FILE *f);

void graph_store(FILE *f, graph *g);

void graph_free(graph *g);

int graph_validate(graph *g);

graph *graph_subgraph(graph *g, int *mask, int *reverse_map);

// Should be called inside parallel region
void graph_subgraph_par(graph *g, graph *sg, int *mask, int *reverse_map, int *forward_map, int *s1, int *s2);