#pragma once
#include <stdio.h>

typedef struct
{
    int N;
    int *V, *E;
    long long *W;
} graph;

graph *graph_parse(FILE *f);

void graph_free(graph *g);

int graph_validate(graph *g);

graph *graph_subgraph(graph *g, int *mask, int *reverse_mapping);