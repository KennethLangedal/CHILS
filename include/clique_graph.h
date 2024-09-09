#pragma once
#include <stdio.h>

typedef struct
{
    int N, NC;   // Number of vertices and cliques
    int *V, *EV; // CSR for vertices to cliques
    int *C, *EC; // CSR for cliques to vertices

    long long *W; // Vertex weights
} clique_graph;

clique_graph *clique_graph_parse_json(FILE *f);

void clique_graph_free(clique_graph *g);

int clique_graph_validate(clique_graph *g);