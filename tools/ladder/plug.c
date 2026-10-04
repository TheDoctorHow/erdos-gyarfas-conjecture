/* plug.c -- hunt frozen 2-plugs: for each cubic H* (graph6 on stdin) and
 * each edge e = pq, let H = H* - e.  Report cases where:
 *   (i)  H has no cycle of length 8 or 16 (n <= 20: covers all == 0 mod 8),
 *   (ii) the set of p-q simple-path residues mod 8 is a singleton {r}.
 * Prints: <g6> edge=p-q r=<r> npaths=<count>
 */
#include "egcore.h"

static EGGraph H;
static int P, Q;
static unsigned char used[EG_MAXV];
static unsigned resmask;
static long long npaths, budget;

static void dfs(int v, int len) {
    if (budget-- < 0) return;
    if (v == Q) { resmask |= 1u << (len % 8); npaths++; return; }
    for (int i = 0; i < H.deg[v]; i++) {
        int w = H.adj[v][i];
        if (used[w]) continue;
        used[w] = 1;
        dfs(w, len + 1);
        used[w] = 0;
        if (budget < 0) return;
    }
}

int main(void) {
    char line[1 << 16];
    long long tot = 0, found = 0;
    while (fgets(line, sizeof line, stdin)) {
        EGGraph G;
        char *nl = strchr(line, '\n'); if (nl) *nl = 0;
        if (!line[0]) continue;
        if (eg_from_graph6(&G, line) != 0) continue;
        tot++;
        for (int u = 0; u < G.n; u++) for (int j = 0; j < G.deg[u]; j++) {
            int v = G.adj[u][j];
            if (v < u) continue;
            /* H = G - uv */
            H = G;
            for (int k = 0; k < H.deg[u]; k++)
                if (H.adj[u][k] == v) { H.adj[u][k] = H.adj[u][--H.deg[u]]; break; }
            for (int k = 0; k < H.deg[v]; k++)
                if (H.adj[v][k] == u) { H.adj[v][k] = H.adj[v][--H.deg[v]]; break; }
            H.is_adj[u][v] = H.is_adj[v][u] = 0;
            /* (i) no 8- or 16-cycles in H */
            if (eg_cycles_of_length(&H, 8, 1, 100000000LL) != 0) continue;
            if (H.n >= 16 && eg_cycles_of_length(&H, 16, 1, 100000000LL) != 0)
                continue;
            /* (ii) frozen p-q residues */
            P = u; Q = v;
            memset(used, 0, sizeof used);
            used[P] = 1;
            resmask = 0; npaths = 0; budget = 200000000LL;
            dfs(P, 0);
            if (budget < 0) { printf("BUDGET %s %d-%d\n", line, u, v); continue; }
            if (npaths > 0 && (resmask & (resmask - 1)) == 0) {
                int r = 0; while (!(resmask & (1u << r))) r++;
                printf("FROZEN %s edge=%d-%d r=%d npaths=%lld\n",
                       line, u, v, r, npaths);
                found++;
            }
        }
    }
    fprintf(stderr, "graphs=%lld frozen-plugs=%lld\n", tot, found);
    return 0;
}
