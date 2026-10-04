/* inflate.c -- triangle-inflation enumeration for the Erdos-Gyarfas hunt.
 *
 * REDUCTION (proved below).  Let G be a connected cubic graph with no C4 and
 * with t >= 1 triangles.  Then:
 *   (i)   the triangles of G are pairwise vertex-disjoint;
 *   (ii)  contracting all of them gives a connected cubic graph H on n-2t
 *         vertices which is TRIANGLE-FREE;
 *   (iii) G is recovered from H by "inflating" a set S of t vertices of H,
 *         each into a triangle;
 *   (iv)  inflate(H,S) is C4-free  <=>  H is triangle-free and every C4 of H
 *         meets S.
 *
 * Proofs.
 * (i)  Two triangles of a cubic graph sharing a vertex x force deg(x) >= 4
 *      unless they share an edge xy; then x ~ y,z,w and y ~ x,z,w, so
 *      z-x-w-y-z is a C4.  Excluded.
 * (ii) The three outside neighbours x',y',z' of a triangle {x,y,z} are
 *      distinct: if x'=y'=w then z-x-w-y-z is a C4.  So the contraction is
 *      cubic and simple.  If H had a triangle through a contracted vertex v,
 *      say v-a-b-v with a~x and b~y (x,y in the triangle), then a-x-y-b-a is
 *      a C4 of G.  Excluded.  A triangle of H avoiding every contracted
 *      vertex would be an uncontracted triangle of G.  Excluded.
 * (iv) A cycle of G meets each triangle it touches in >= 2 vertices (a cycle
 *      through x uses two of x's three edges, two of which stay inside).
 *      Enumerating the cases, the only C4s of G are the C4s of H that avoid S.
 *
 * Consequence: enumerating cubic C4-free graphs on n vertices reduces to
 *   t = 0 : cubic graphs of girth >= 5 on n vertices, plus
 *   t >= 1: triangle-free cubic graphs on n-2t vertices with a t-element
 *           hitting set for their C4s.
 * The bases are far scarcer than cubic graphs on n vertices, and this is
 * exactly where the known near-counterexamples live (Markstrom's four
 * 24-vertex graphs have 3, 4, 6 and 7 triangles, i.e. bases on 18, 16, 12
 * and 10 vertices).
 *
 * Usage:  snarkhunter <m> 4 o g s | inflate <t> [-v]
 *         reads triangle-free cubic graphs H on m vertices, emits statistics
 *         and any counterexample on n = m + 2t vertices.
 */
#include "egcore.h"

static int T;                 /* number of vertices to inflate */
static long long n_base = 0, n_cand = 0, n_noC4 = 0, n_noC4C8 = 0, n_ctr = 0;
static long long histo[8];
static FILE *nearf = NULL;

static int c4v[512][4];       /* the C4s of H, as vertex quadruples */
static int nc4;

static void find_C4s(const EGGraph *H) {
    int u, i, j, x, m = H->n;
    nc4 = 0;
    for (u = 0; u < m; u++)
        for (i = 0; i < H->deg[u]; i++)
            for (j = i + 1; j < H->deg[u]; j++) {
                int a = H->adj[u][i], b = H->adj[u][j];
                for (x = 0; x < H->deg[a]; x++) {
                    int w = H->adj[a][x];
                    if (w == u || !H->is_adj[w][b]) continue;
                    if (w < u) continue;           /* canonical: u smallest-ish */
                    if (nc4 < 512) {
                        c4v[nc4][0] = u; c4v[nc4][1] = a; c4v[nc4][2] = w; c4v[nc4][3] = b;
                        nc4++;
                    }
                }
            }
}

/* Build inflate(H,S) into G.  S given as a bitmask over H's vertices. */
static void build_inflation(const EGGraph *H, unsigned long long S, EGGraph *G) {
    int m = H->n, i, p, idx = 0;
    int base[EG_MAXV];        /* first new vertex index for an inflated vertex */
    int img[EG_MAXV];         /* image of a non-inflated vertex */
    for (i = 0; i < m; i++) {
        if ((S >> i) & 1ULL) { base[i] = idx; idx += 3; img[i] = -1; }
        else { img[i] = idx++; base[i] = -1; }
    }
    eg_clear(G, idx);
    /* triangles */
    for (i = 0; i < m; i++)
        if (base[i] >= 0) {
            eg_add_edge(G, base[i] + 0, base[i] + 1);
            eg_add_edge(G, base[i] + 1, base[i] + 2);
            eg_add_edge(G, base[i] + 2, base[i] + 0);
        }
    /* original edges: edge (i, adj[i][p]) attaches at slot p of i's triangle */
    for (i = 0; i < m; i++)
        for (p = 0; p < H->deg[i]; p++) {
            int j = H->adj[i][p], q;
            if (j < i) continue;
            for (q = 0; q < H->deg[j]; q++) if (H->adj[j][q] == i) break;
            {
                int a = (base[i] >= 0) ? base[i] + p : img[i];
                int b = (base[j] >= 0) ? base[j] + q : img[j];
                eg_add_edge(G, a, b);
            }
        }
}

static void test_candidate(const EGGraph *H, unsigned long long S) {
    EGGraph G;
    unsigned mask, unk;
    int k, lowest = -1;
    n_cand++;
    build_inflation(H, S, &G);
    if (eg_has_C4(&G)) return;                /* safety net: should never fire */
    n_noC4++;
    if (eg_cycles_of_length(&G, 8, 1, 100000000LL) > 0) return;
    n_noC4C8++;
    mask = eg_pow2_mask(&G, 400000000LL, &unk);
    for (k = 2; k < 32; k++) if (mask & (1u << k)) { lowest = k; break; }
    if (lowest >= 0 && lowest < 8) histo[lowest]++;
    {
        char g6[512];
        eg_to_graph6(&G, g6);
        if (nearf) { fprintf(nearf, "%s n=%d mask=%u unk=%u lowest=%d\n", g6, G.n, mask, unk,
                             lowest >= 0 ? (1 << lowest) : -1); fflush(nearf); }
        if (mask == 0 && unk == 0) {
            n_ctr++;
            printf("COUNTEREXAMPLE %s n=%d\n", g6, G.n); fflush(stdout);
            fprintf(stderr, "\n*** COUNTEREXAMPLE *** %s n=%d\n", g6, G.n);
        }
    }
}

/* Enumerate hitting sets S of size exactly T for the C4s of H, then test. */
static void hit(const EGGraph *H, unsigned long long S, int chosen, int firstfree) {
    int i, c;
    /* find an unhit C4 */
    int unhit = -1;
    for (c = 0; c < nc4; c++) {
        unsigned long long q = (1ULL << c4v[c][0]) | (1ULL << c4v[c][1]) |
                               (1ULL << c4v[c][2]) | (1ULL << c4v[c][3]);
        if (!(S & q)) { unhit = c; break; }
    }
    if (unhit >= 0) {
        if (chosen == T) return;                       /* cannot hit it any more */
        for (i = 0; i < 4; i++) {
            int v = c4v[unhit][i];
            if ((S >> v) & 1ULL) continue;
            if (v < firstfree) continue;               /* avoid duplicate sets */
            hit(H, S | (1ULL << v), chosen + 1, 0);
        }
        return;
    }
    /* all C4s hit: pad up to size T with any remaining vertices */
    if (chosen == T) { test_candidate(H, S); return; }
    for (i = firstfree; i < H->n; i++) {
        if ((S >> i) & 1ULL) continue;
        hit(H, S | (1ULL << i), chosen + 1, i + 1);
    }
}

int main(int argc, char **argv) {
    char line[1 << 16];
    const char *nf = getenv("EG_NEARFILE");
    int i;
    T = (argc > 1) ? atoi(argv[1]) : 1;
    if (nf) nearf = fopen(nf, "a");
    while (fgets(line, sizeof(line), stdin)) {
        EGGraph H;
        char *p = strchr(line, '\n'); if (p) *p = 0;
        p = strchr(line, '\r'); if (p) *p = 0;
        if (!line[0]) continue;
        if (eg_from_graph6(&H, line) != 0) { fprintf(stderr, "BAD %s\n", line); continue; }
        n_base++;
        find_C4s(&H);
        hit(&H, 0ULL, 0, 0);
    }
    fprintf(stderr, "INFSTAT t=%d bases=%lld candidates=%lld noC4=%lld noC4C8=%lld ctr=%lld\n",
            T, n_base, n_cand, n_noC4, n_noC4C8, n_ctr);
    for (i = 0; i < 8; i++)
        if (histo[i]) fprintf(stderr, "INFSTAT   among noC4C8: smallest 2^k = %d : %lld\n", 1 << i, histo[i]);
    if (nearf) fclose(nearf);
    return 0;
}
