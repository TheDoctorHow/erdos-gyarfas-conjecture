/* core4.c -- k=4 tube-frame residue scan (straight twist, pi = id).
 * Variables mod 8: gaps g[4] on C1 (sum==4), h[4] on C2 (sum==4),
 * path residues p[4].  Webs:
 *  - pair-webs (i,j): p_i+p_j+A+B, A in {a, 4-a} (a = gap-sum u_i..u_j),
 *    B likewise on C2; all four combos must be != 0 (mod 8).
 *  - 4-webs (two patterns): p-sum + {g1+g3+h2+h4, g2+g4+h1+h3} != 0.
 * Avoider demand: for each C1-gap i: some web == 4 whose C1-arcs avoid
 * gap i; same for C2-gaps.  Count survivors.
 */
#include <stdio.h>

static int gsum(const int *g, int i, int j) {   /* gaps i..j-1 cyclic */
    int s = 0;
    for (int k = i; k != j; k = (k + 1) & 3) s += g[k];
    return s & 7;
}

int main(void) {
    long long total = 0, survive = 0, av_survive = 0;
    int g[4], h[4], p[4];
    for (g[0] = 0; g[0] < 8; g[0]++)
    for (g[1] = 0; g[1] < 8; g[1]++)
    for (g[2] = 0; g[2] < 8; g[2]++) {
        g[3] = (4 - g[0] - g[1] - g[2]) & 7;
        for (h[0] = 0; h[0] < 8; h[0]++)
        for (h[1] = 0; h[1] < 8; h[1]++)
        for (h[2] = 0; h[2] < 8; h[2]++) {
            h[3] = (4 - h[0] - h[1] - h[2]) & 7;
            for (p[0] = 0; p[0] < 8; p[0]++)
            for (p[1] = 0; p[1] < 8; p[1]++)
            for (p[2] = 0; p[2] < 8; p[2]++)
            for (p[3] = 0; p[3] < 8; p[3]++) {
                total++;
                int ok = 1;
                int has4[8];            /* webs == 4 avoiding C1-gap i (bit i), C2-gap i (bit 4+i) */
                for (int i = 0; i < 8; i++) has4[i] = 0;
                /* pair-webs */
                for (int i = 0; i < 4 && ok; i++)
                for (int j = i + 1; j < 4 && ok; j++) {
                    int a = gsum(g, i, j), b = gsum(h, i, j);
                    int A[2] = {a, (4 - a) & 7};
                    int B[2] = {b, (4 - b) & 7};
                    for (int x = 0; x < 2 && ok; x++)
                    for (int y = 0; y < 2 && ok; y++) {
                        int L = (p[i] + p[j] + A[x] + B[y]) & 7;
                        if (L == 0) { ok = 0; break; }
                        if (L == 4) {
                            /* which C1 gaps does this web avoid?  arc x=0 uses
                               gaps i..j-1; x=1 uses the complement */
                            for (int q = 0; q < 4; q++) {
                                int inarc = 0;
                                for (int k = i; k != j; k = (k + 1) & 3)
                                    if (k == q) inarc = 1;
                                if ((x == 0) ? !inarc : inarc) has4[q] |= 1;
                            }
                            for (int q = 0; q < 4; q++) {
                                int inarc = 0;
                                for (int k = i; k != j; k = (k + 1) & 3)
                                    if (k == q) inarc = 1;
                                if ((y == 0) ? !inarc : inarc) has4[4 + q] |= 1;
                            }
                        }
                    }
                }
                if (!ok) continue;
                /* 4-webs (straight) */
                int ps = (p[0] + p[1] + p[2] + p[3]) & 7;
                int w1 = (ps + g[0] + g[2] + h[1] + h[3]) & 7;
                int w2 = (ps + g[1] + g[3] + h[0] + h[2]) & 7;
                if (w1 == 0 || w2 == 0) continue;
                survive++;
                /* avoider demand */
                int av = 1;
                for (int q = 0; q < 8 && av; q++)
                    if (!has4[q]) av = 0;
                if (av) av_survive++;
            }
        }
    }
    printf("total=%lld  web-safe=%lld  web-safe+avoiders=%lld\n",
           total, survive, av_survive);
    return 0;
}
