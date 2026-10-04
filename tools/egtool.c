/* egtool.c -- read graph6 from stdin, report power-of-two cycle structure.
 *
 *   egtool           : one line per graph:  <g6> n=<n> pow2=<list> min2k=<v> [FLAG]
 *   egtool -c        : also print exact counts of each 2^k cycle
 *   egtool -s        : print the full cycle spectrum (slow, small graphs)
 *   egtool -q        : quiet; only print COUNTEREXAMPLE / NEARMISS lines
 *
 * A COUNTEREXAMPLE line means: min degree >= 3 and no cycle of length 2^k.
 * A NEARMISS line means: no C4 and no C8 (but some larger 2^k present).
 */
#include "egcore.h"

int main(int argc, char **argv) {
    char line[1 << 16];
    int want_counts = 0, want_spec = 0, quiet = 0, i;
    long long budget = 400000000LL;
    long long ngraphs = 0, nnear = 0, nctr = 0, nunknown = 0;
    long long histo[32];

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-c")) want_counts = 1;
        else if (!strcmp(argv[i], "-s")) want_spec = 1;
        else if (!strcmp(argv[i], "-q")) quiet = 1;
        else if (!strcmp(argv[i], "-b")) budget = atoll(argv[++i]);
    }
    memset(histo, 0, sizeof(histo));

    while (fgets(line, sizeof(line), stdin)) {
        EGGraph G;
        unsigned mask, unk;
        int mind, k, rc;
        char *nl = strchr(line, '\n'); if (nl) *nl = 0;
        nl = strchr(line, '\r'); if (nl) *nl = 0;
        if (line[0] == 0) continue;
        rc = eg_from_graph6(&G, line);
        if (rc != 0) { fprintf(stderr, "BAD graph6 (rc=%d): %s\n", rc, line); continue; }
        ngraphs++;

        mind = 1 << 20;
        for (k = 0; k < G.n; k++) if (G.deg[k] < mind) mind = G.deg[k];

        mask = eg_pow2_mask(&G, budget, &unk);
        if (unk) nunknown++;

        {
            int lowest = -1;
            for (k = 2; k < 32; k++) if (mask & (1u << k)) { lowest = k; break; }
            if (lowest >= 0 && lowest < 32) histo[lowest]++;

            if (mind >= 3 && mask == 0 && unk == 0) {
                nctr++;
                printf("COUNTEREXAMPLE %s n=%d mindeg=%d\n", line, G.n, mind);
                fflush(stdout);
            } else if (mind >= 3 && unk) {
                printf("UNKNOWN %s n=%d unkmask=%u knownmask=%u\n", line, G.n, unk, mask);
                fflush(stdout);
            } else if (!(mask & (1u << 2)) && !(mask & (1u << 3))) {
                nnear++;
                printf("NEARMISS %s n=%d mindeg=%d mask=%u lowest2k=%d\n",
                       line, G.n, mind, mask, lowest >= 0 ? (1 << lowest) : -1);
                fflush(stdout);
            }
        }

        if (!quiet) {
            printf("%s n=%d mindeg=%d pow2=", line, G.n, mind);
            for (k = 2; (1 << k) <= G.n; k++) if (mask & (1u << k)) printf("%d,", 1 << k);
            if (want_counts) {
                printf(" counts=");
                for (k = 2; (1 << k) <= G.n; k++)
                    printf("%d:%lld,", 1 << k, eg_cycles_of_length(&G, 1 << k, 0, budget));
            }
            if (want_spec) {
                printf(" spec=");
                for (k = 3; k <= G.n; k++)
                    if (eg_cycles_of_length(&G, k, 1, budget) > 0) printf("%d,", k);
            }
            printf("\n");
        }
    }
    fprintf(stderr, "graphs=%lld nearmiss=%lld counterexamples=%lld unknown=%lld\n",
            ngraphs, nnear, nctr, nunknown);
    fprintf(stderr, "histogram of smallest power-of-two cycle:\n");
    for (i = 0; i < 32; i++) if (histo[i]) fprintf(stderr, "  2^%d=%d : %lld\n", i, 1 << i, histo[i]);
    return 0;
}
