#pragma once
#include <stdio.h>

typedef struct
{
    int n;
    int *V, *E;
    long long *W;
} graph;

graph *graph_parse(FILE *f);

void graph_store(FILE *f, graph *g);

void graph_free(graph *g);

int graph_validate(graph *g);

graph *graph_subgraph(graph *g, int *Mask, int *RM);

// Should be called inside parallel region
void graph_subgraph_par(graph *g, graph *sg, int *Mask, int *RM, int *FM, int *S1, int *S2);