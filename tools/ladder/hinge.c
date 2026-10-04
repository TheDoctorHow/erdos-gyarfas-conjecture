/* hinge.c -- does the plug operad cover the graphs whose spectrum has holes?
 *
 * A 3-terminal plug substituted at a vertex of a cubic graph creates a NON-TRIVIAL
 * 3-edge-cut (both sides having >= 2 vertices), and conversely every non-trivial
 * 3-edge-cut of a cubic graph exhibits one side as a plug.  For a 3-connected cubic
 * graph that is the same as saying its cyclic edge connectivity equals 3.
 *
 * So "G decomposes through the plug operad"  <=>  G has a non-trivial 3-edge-cut.
 *
 * This program cross-tabulates, over 3-connected cubic graphs:
 *      hole-free spectrum?   x   decomposable (cyc-3)?   x   dyadically fat?
 * The class that matters is HOLED and NON-DECOMPOSABLE: those are the graphs the
 * operad cannot reach, and they are what Gap 3 would still have to handle.
 */
#include "egcore.h"
#include <stdint.h>

static int ea[64], eb[64], M;
static int par[EG_MAXV];
static int find(int x){ while(par[x]!=x){ par[x]=par[par[x]]; x=par[x]; } return x; }

/* number of vertices reachable from 0 with the listed edge indices removed */
static int reach(const EGGraph *G, int x, int y, int z) {
    int seen[EG_MAXV] = {0}, st[EG_MAXV], sp = 0, cnt = 1;
    seen[0] = 1; st[sp++] = 0;
    while (sp) {
        int u = st[--sp];
        for (int j = 0; j < G->deg[u]; j++) {
            int v = G->adj[u][j];
            int skip = 0;
            for (int t = 0; t < 3; t++) {
                int e = (t == 0) ? x : (t == 1) ? y : z;
                if (e < 0) continue;
                if ((ea[e] == u && eb[e] == v) || (eb[e] == u && ea[e] == v)) skip = 1;
            }
            if (skip || seen[v]) continue;
            seen[v] = 1; st[sp++] = v; cnt++;
        }
    }
    return cnt;
}

int main(void) {
    char line[1 << 18];
    long long tot = 0, elig = 0;
    long long a[2][2][2] = {{{0}}};   /* [holefree][decomposable][fat] */
    int shown = 0;
    while (fgets(line, sizeof line, stdin)) {
        EGGraph G;
        char *p = strchr(line, '\n'); if (p) *p = 0;
        if (!line[0]) continue;
        if (eg_from_graph6(&G, line) != 0) continue;
        int n = G.n, mind = 99;
        for (int i = 0; i < n; i++) if (G.deg[i] < mind) mind = G.deg[i];
        if (mind != 3) continue;
        tot++;
        M = 0;
        for (int u = 0; u < n; u++) for (int j = 0; j < G.deg[u]; j++) {
            int v = G.adj[u][j]; if (v > u) { ea[M] = u; eb[M] = v; M++; }
        }
        /* 3-edge-connected (= 3-connected for cubic) */
        int ok3 = 1;
        for (int x = 0; x < M && ok3; x++) for (int y = x; y < M && ok3; y++)
            if (reach(&G, x, y, -1) < n) ok3 = 0;
        if (!ok3) continue;
        elig++;

        /* non-trivial 3-edge-cut?  (both sides >= 2 vertices) */
        int decomp = 0;
        for (int x = 0; x < M && !decomp; x++)
            for (int y = x + 1; y < M && !decomp; y++)
                for (int z = y + 1; z < M && !decomp; z++) {
                    int r = reach(&G, x, y, z);
                    if (r < n && r >= 2 && n - r >= 2) decomp = 1;
                }

        /* spectrum via the GF(2) cycle space (= all 2-regular subgraphs, cubic case) */
        int inT[64] = {0}, pare[EG_MAXV], parv[EG_MAXV], seen[EG_MAXV] = {0}, order[EG_MAXV], no = 0;
        for (int i = 0; i < n; i++) { pare[i] = -1; parv[i] = -1; }
        seen[0] = 1; order[no++] = 0;
        for (int qi = 0; qi < no; qi++) {
            int u = order[qi];
            for (int j = 0; j < G.deg[u]; j++) {
                int v = G.adj[u][j]; if (seen[v]) continue;
                for (int t = 0; t < M; t++)
                    if ((ea[t] == u && eb[t] == v) || (ea[t] == v && eb[t] == u)) { inT[t] = 1; pare[v] = t; break; }
                parv[v] = u; seen[v] = 1; order[no++] = v;
            }
        }
        int depth[EG_MAXV]; depth[0] = 0;
        for (int qi = 1; qi < no; qi++) { int v = order[qi]; depth[v] = depth[parv[v]] + 1; }
        uint64_t fc[40]; int nf = 0;
        for (int t = 0; t < M; t++) {
            if (inT[t]) continue;
            uint64_t msk = 1ULL << t; int u = ea[t], v = eb[t];
            while (u != v) {
                if (depth[u] < depth[v]) { msk ^= 1ULL << pare[v]; v = parv[v]; }
                else                     { msk ^= 1ULL << pare[u]; u = parv[u]; }
            }
            fc[nf++] = msk;
        }
        static char isCyc[64];
        for (int i = 0; i <= M; i++) isCyc[i] = 0;
        uint64_t cur = 0; long long total = 1LL << nf;
        for (long long g = 1; g < total; g++) {
            cur ^= fc[__builtin_ctzll((uint64_t)g)];
            int sz = __builtin_popcountll(cur);
            if (!sz || isCyc[sz]) continue;
            for (int i = 0; i < n; i++) par[i] = i;
            uint64_t w = cur;
            while (w) { int t = __builtin_ctzll(w); w &= w - 1; int q = find(ea[t]), r = find(eb[t]); if (q != r) par[q] = r; }
            int mark[EG_MAXV] = {0}, comp = 0; w = cur;
            while (w) { int t = __builtin_ctzll(w); w &= w - 1; mark[find(ea[t])] = 1; }
            for (int i = 0; i < n; i++) if (mark[i]) comp++;
            if (comp == 1) isCyc[sz] = 1;
        }

        int gg = 0, cc = 0, bip = 1;
        for (int L = 3; L <= M; L++) if (isCyc[L]) { if (!gg) gg = L; cc = L; if (L & 1) bip = 0; }
        int holefree = 1, step = bip ? 2 : 1;
        for (int L = gg; L <= cc; L += step) if (!isCyc[L]) holefree = 0;

        int fat = 0;
        for (int d = 1; d <= 2; d++)
            for (int a0 = 3; a0 <= M; a0++) {
                if (!isCyc[a0]) continue;
                if (d == 2 && (a0 & 1)) continue;
                int b = a0; while (b + d <= M && isCyc[b + d]) b += d;
                if (b >= 2 * a0 - d) fat = 1;
            }

        a[holefree][decomp][fat]++;
        if (!holefree && !decomp && shown < 12) {
            shown++;
            printf("HOLED + NOT DECOMPOSABLE%s : %s n=%d g=%d c=%d spec=",
                   fat ? "" : "  *** AND NOT FAT ***", line, n, gg, cc);
            for (int L = 3; L <= M; L++) if (isCyc[L]) printf("%d,", L);
            printf("\n");
        }
    }
    printf("cubic=%lld  3-connected=%lld\n", tot, elig);
    printf("                         decomposable(cyc-3)      NOT decomposable(cyc>=4)\n");
    printf("  hole-free  fat        %12lld              %12lld\n", a[1][1][1], a[1][0][1]);
    printf("  hole-free  NOT fat    %12lld              %12lld\n", a[1][1][0], a[1][0][0]);
    printf("  HOLED      fat        %12lld              %12lld\n", a[0][1][1], a[0][0][1]);
    printf("  HOLED      NOT fat    %12lld              %12lld\n", a[0][1][0], a[0][0][0]);
    return 0;
}
