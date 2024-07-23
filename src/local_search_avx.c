#include "local_search_avx.h"

#include <stdlib.h>
#include <immintrin.h>

local_search_avx local_search_avx_init(graph g)
{
    local_search_avx ls;
    ls.N = g.N;
    ls.C = malloc(sizeof(int));
    *ls.C = 0;

    ls.IS = aligned_alloc(32, sizeof(int) * (g.N + 1));
    for (int i = 0; i < g.N + 1; i++)
        ls.IS[i] = 0;

    ls.V = malloc(sizeof(int) * (g.N + 1));
    ls.V[0];
    for (int i = 0; i < g.N; i++)
        ls.V[i + 1] = ls.V[i] + (((g.V[i + 1] - g.V[i]) + 7) & (~7));

    ls.E = aligned_alloc(32, sizeof(int) * ls.V[g.N]);
    ls.W = aligned_alloc(32, sizeof(int) * (g.N + 1));
    for (int i = 0; i < g.N; i++)
    {
        ls.W[i] = g.W[i];
        int j = g.V[i], k = ls.V[i];
        while (j < g.V[i + 1])
            ls.E[k++] = g.E[j++];
        while (k < ls.V[i + 1])
            ls.E[k++] = g.N;
    }
    ls.W[g.N] = 0;

    return ls;
}

void local_search_avx_free(local_search_avx ls)
{
    free(ls.C);
    free(ls.IS);
    free(ls.V);
    free(ls.E);
    free(ls.W);
}

void local_search_avx_add_vertex(local_search_avx ls, int u)
{
    if (ls.IS[u])
        return;

    ls.IS[u] = 0xffffffffu;
    *ls.C += ls.W[u];

    for (int i = ls.V[u]; i < ls.V[u + 1]; i++)
    {
        int v = ls.E[i];
        if (ls.IS[v])
        {
            ls.IS[v] = 0;
            *ls.C -= ls.W[v];
        }
    }
}

void local_search_avx_remove_vertex(local_search_avx ls, int u)
{
    if (!ls.IS[u])
        return;

    ls.IS[u] = 0;
    *ls.C -= ls.W[u];
}

void local_search_avx_greedy(local_search_avx ls, int *order)
{
    int imp = 1;
    while (imp)
    {
        imp = 0;
        for (int i = 0; i < ls.N; i++)
        {
            int u = order[i];
            if (ls.IS[u])
                continue;

            __m256i v_nw = _mm256_setzero_si256();
            for (int j = ls.V[u]; j < ls.V[u + 1]; j += 8)
            {
                __m256i index = _mm256_load_si256((__m256i const *)(ls.E + j));
                __m256i v_is = _mm256_i32gather_epi32(ls.IS, index, 4);
                __m256i v_w = _mm256_i32gather_epi32(ls.W, index, 4);

                v_nw = _mm256_add_epi32(v_nw, _mm256_and_si256(v_is, v_w));
            }

            v_nw = _mm256_hadd_epi32(v_nw, v_nw);
            v_nw = _mm256_hadd_epi32(v_nw, v_nw);

            int nw = _mm256_extract_epi32(v_nw, 0) + _mm256_extract_epi32(v_nw, 4);

            if (nw < ls.W[u])
            {
                local_search_avx_add_vertex(ls, u);
                imp = 1;
            }
        }
    }
}