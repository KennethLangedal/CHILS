#pragma once
#include <stdio.h>

typedef struct
{
    int N;
    int *V, *E;
    int *W;
} graph;

graph graph_parse(FILE *f);

void graph_free(graph g);

int graph_validate(int N, const int *V, const int *E);