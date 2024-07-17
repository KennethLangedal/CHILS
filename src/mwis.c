#include "mwis.h"

int mwis_validate(graph g, int *IS)
{
    int C = 0;
    for (int u = 0; u < g.N; u++)
    {
        if (!IS[u])
            continue;

        C += g.W[u];
        for (int i = g.V[u]; i < g.V[u + 1]; i++)
        {
            int v = g.E[i];
            if (IS[v])
                return -1;
        }
    }
    return C;
}