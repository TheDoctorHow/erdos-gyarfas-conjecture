/* egcore.h -- core machinery for the Erdos-Gyarfas hunt.
 *
 * Exact detection/counting of cycles of a GIVEN length in a sparse graph.
 * Written to be independently auditable: no clever tricks, one algorithm,
 * heavy assertions available under -DEG_DEBUG.
 *
 * Method for "is there a cycle of length exactly L?":
 *   Every cycle has a unique minimum vertex r.  For each r we search for
 *   simple paths r = p_0, p_1, ..., p_{L-1} with all p_i > r (i>0),
 *   closed by the edge p_{L-1}p_0.  A BFS distance array from r inside the
 *   induced subgraph on {r, r+1, ..., n-1} prunes: from position k we still
 *   need to reach r within L-k edges, so prune when dist[cur] > L-k.
 *   To avoid counting each cycle twice we require p_1 < p_{L-1} when counting.
 *
 * This is exact (no heuristics) but exponential in L; a node budget lets the
 * caller detect the rare cases where it would be too slow.
 */
#ifndef EGCORE_H
#define EGCORE_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef EG_MAXV
#define EG_MAXV 128
#endif
#ifndef EG_MAXD
#define EG_MAXD 16
#endif

typedef struct {
    int n;
    int deg[EG_MAXV];
    int adj[EG_MAXV][EG_MAXD];
    unsigned char is_adj[EG_MAXV][EG_MAXV];
} EGGraph;

static void eg_clear(EGGraph *G, int n) {
    int i;
    G->n = n;
    for (i = 0; i < n; i++) { G->deg[i] = 0; memset(G->is_adj[i], 0, (size_t) n); }
}

/* returns 0 on success, -1 if the edge is a loop/duplicate/overflow */
static int eg_add_edge(EGGraph *G, int u, int v) {
    if (u == v) return -1;
    if (G->is_adj[u][v]) return -1;
    if (G->deg[u] >= EG_MAXD || G->deg[v] >= EG_MAXD) return -1;
    G->adj[u][G->deg[u]++] = v;
    G->adj[v][G->deg[v]++] = u;
    G->is_adj[u][v] = G->is_adj[v][u] = 1;
    return 0;
}

/* ---- search state (single-threaded per translation unit) ---- */
typedef struct {
    const EGGraph *G;
    int L;              /* target cycle length */
    int root;
    int dist[EG_MAXV];  /* BFS dist from root within {root..n-1} */
    unsigned char used[EG_MAXV];
    int path[EG_MAXV];
    long long budget;   /* decremented per DFS node; <0 means exhausted */
    long long count;    /* cycles found (when counting) */
    int stop_at_first;
    int found;
} EGSearch;

static void eg_bfs_from_root(EGSearch *S) {
    const EGGraph *G = S->G;
    int q[EG_MAXV], head = 0, tail = 0, i, r = S->root;
    for (i = 0; i < G->n; i++) S->dist[i] = 1 << 20;
    S->dist[r] = 0;
    q[tail++] = r;
    while (head < tail) {
        int u = q[head++], j;
        for (j = 0; j < G->deg[u]; j++) {
            int w = G->adj[u][j];
            if (w < r) continue;                 /* restrict to {root..n-1} */
            if (S->dist[w] > S->dist[u] + 1) {
                S->dist[w] = S->dist[u] + 1;
                q[tail++] = w;
            }
        }
    }
}

/* cur is at index k in the path (path[0]==root). */
static void eg_dfs(EGSearch *S, int cur, int k) {
    const EGGraph *G = S->G;
    int j;
    if (S->budget-- < 0) return;
    if (S->found && S->stop_at_first) return;

    if (k == S->L - 1) {
        if (G->is_adj[cur][S->root]) {
            /* orientation tie-break so each cycle is counted once */
            if (S->L == 3 || S->path[1] < cur) { S->count++; S->found = 1; }
        }
        return;
    }
    /* need to get back to root in L-k edges */
    if (S->dist[cur] > S->L - k) return;

    for (j = 0; j < G->deg[cur]; j++) {
        int w = G->adj[cur][j];
        if (w <= S->root) continue;
        if (S->used[w]) continue;
        /* cannot close early: if w is adjacent to root we may still continue,
           the length check at k==L-1 handles correctness */
        S->used[w] = 1;
        S->path[k + 1] = w;
        eg_dfs(S, w, k + 1);
        S->used[w] = 0;
        if (S->found && S->stop_at_first) return;
        if (S->budget < 0) return;
    }
}

/* Count (or detect) cycles of length exactly L.
 * Returns the count, or -1 if the node budget was exhausted (answer unknown).
 * If stop_at_first, the return value is 0 or 1 (1 = exists).            */
static long long eg_cycles_of_length(const EGGraph *G, int L, int stop_at_first,
                                     long long budget) {
    EGSearch S;
    int r;
    if (L < 3 || L > G->n) return 0;
    S.G = G; S.L = L; S.budget = budget; S.count = 0;
    S.stop_at_first = stop_at_first; S.found = 0;
    for (r = 0; r + L <= G->n; r++) {
        if (G->deg[r] < 2) continue;
        S.root = r;
        eg_bfs_from_root(&S);
        memset(S.used, 0, sizeof(S.used));
        S.used[r] = 1;
        S.path[0] = r;
        S.found = 0;
        eg_dfs(&S, r, 0);
        if (S.budget < 0) return -1;
        if (stop_at_first && S.count > 0) return 1;
    }
    return S.count;
}

/* ---------- powers of two ---------- */

/* Returns a bitmask: bit k set  <=>  G has a cycle of length 2^k.
 * Only k with 4 <= 2^k <= n are examined (2^0,2^1 impossible; 2^k>=4).
 * If any test runs out of budget, *unknown is set to the mask of
 * undecided k and those bits are NOT set in the result.               */
static unsigned eg_pow2_mask(const EGGraph *G, long long budget, unsigned *unknown) {
    unsigned mask = 0, unk = 0;
    int k;
    for (k = 2; (1 << k) <= G->n; k++) {
        long long r = eg_cycles_of_length(G, 1 << k, 1, budget);
        if (r < 0) unk |= 1u << k;
        else if (r > 0) mask |= 1u << k;
    }
    if (unknown) *unknown = unk;
    return mask;
}

/* Fast, allocation-free C4 test: two vertices with >=2 common neighbours. */
static int eg_has_C4(const EGGraph *G) {
    int u, i, j, w, x;
    for (u = 0; u < G->n; u++) {
        for (i = 0; i < G->deg[u]; i++) {
            for (j = i + 1; j < G->deg[u]; j++) {
                int a = G->adj[u][i], b = G->adj[u][j];
                /* common neighbour of a,b other than u  => C4 */
                for (x = 0; x < G->deg[a]; x++) {
                    w = G->adj[a][x];
                    if (w != u && G->is_adj[w][b]) return 1;
                }
            }
        }
    }
    return 0;
}

/* ---------- graph6 I/O ---------- */

/* Decode a graph6 line into G. Returns 0 on success. */
static int eg_from_graph6(EGGraph *G, const char *s) {
    int n, i, j, k, bitpos;
    const unsigned char *p = (const unsigned char *)s;
    if (*p == '>') {                     /* skip optional >>graph6<< header */
        while (*p && *p != '\n' && *p != ':' && *p != '?' ) p++;
    }
    if (*p < 63) return -1;
    if (*p == 126) {                     /* n >= 63 */
        p++;
        if (*p == 126) {                 /* huge */
            p++;
            n = 0;
            for (i = 0; i < 6; i++) { n = (n << 6) | (*p++ - 63); }
        } else {
            n = 0;
            for (i = 0; i < 3; i++) { n = (n << 6) | (*p++ - 63); }
        }
    } else {
        n = *p++ - 63;
    }
    if (n > EG_MAXV) return -2;
    eg_clear(G, n);
    bitpos = 0;
    for (j = 1; j < n; j++) {
        for (i = 0; i < j; i++) {
            int byte = bitpos / 6, bit = 5 - (bitpos % 6);
            int c = p[byte];
            if (c < 63) return -3;
            if (((c - 63) >> bit) & 1) {
                if (eg_add_edge(G, i, j) != 0) return -4;
            }
            bitpos++;
        }
    }
    (void)k;
    return 0;
}

static void eg_to_graph6(const EGGraph *G, char *out) {
    int n = G->n, i, j, bitpos = 0, len, t;
    unsigned char bits[EG_MAXV * EG_MAXV / 2 / 6 + 8];
    char *o = out;
    memset(bits, 0, sizeof(bits));
    for (j = 1; j < n; j++)
        for (i = 0; i < j; i++) {
            if (G->is_adj[i][j]) bits[bitpos / 6] |= 1 << (5 - (bitpos % 6));
            bitpos++;
        }
    if (n <= 62) {
        *o++ = (char)(n + 63);
    } else {
        *o++ = 126;
        *o++ = (char)(((n >> 12) & 63) + 63);
        *o++ = (char)(((n >> 6) & 63) + 63);
        *o++ = (char)((n & 63) + 63);
    }
    len = (bitpos + 5) / 6;
    for (t = 0; t < len; t++) *o++ = (char)(bits[t] + 63);
    *o = 0;
}

#endif /* EGCORE_H */
