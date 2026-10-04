/* g5filter.c -- snarkhunter plugin: Erdos-Gyarfas counterexample hunt. */
#include "nausparse.h"
#include "snarkhunter.h"

#define FILTER  eg_filter
#define SUMMARY eg_summary
/* detail:
 *
 * Called on every completed cubic graph.  Rejects (returns 0) any graph
 * that has a cycle of length a power of two -- i.e. it keeps ONLY
 * counterexamples.  Graphs with no C4 and no C8 are logged as near-misses.
 *
 * Env vars:
 *   EG_NEARFILE   path for the near-miss catalogue (default eg_nearmiss.txt)
 *   EG_TAG        string prefixed to summary lines
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define EGF_MAXV 64

static unsigned long long egf_seen = 0;
static unsigned long long egf_noC4 = 0;
static unsigned long long egf_noC4C8 = 0;
static unsigned long long egf_ctr = 0;
static unsigned long long egf_histo[8];      /* smallest 2^k present, index k */
static FILE *egf_near = NULL;

/* --- local cubic-graph cycle machinery (degree exactly 3) --- */
static int   egf_n;
static unsigned char egf_adj[EGF_MAXV][3];
static unsigned long long egf_nbrmask[EGF_MAXV];
static int   egf_dist[EGF_MAXV];
static unsigned long long egf_used;
static int   egf_L, egf_root, egf_found;

static void egf_bfs(int r) {
    int q[EGF_MAXV], head = 0, tail = 0, i;
    for (i = 0; i < egf_n; i++) egf_dist[i] = 1000;
    egf_dist[r] = 0; q[tail++] = r;
    while (head < tail) {
        int u = q[head++], j;
        for (j = 0; j < 3; j++) {
            int w = egf_adj[u][j];
            if (w < r) continue;
            if (egf_dist[w] > egf_dist[u] + 1) { egf_dist[w] = egf_dist[u] + 1; q[tail++] = w; }
        }
    }
}

static void egf_dfs(int cur, int k) {
    int j;
    if (egf_found) return;
    if (k == egf_L - 1) {
        if ((egf_nbrmask[cur] >> egf_root) & 1ULL) egf_found = 1;
        return;
    }
    if (egf_dist[cur] > egf_L - k) return;
    for (j = 0; j < 3; j++) {
        int w = egf_adj[cur][j];
        if (w <= egf_root) continue;
        if ((egf_used >> w) & 1ULL) continue;
        egf_used |= 1ULL << w;
        egf_dfs(w, k + 1);
        egf_used &= ~(1ULL << w);
        if (egf_found) return;
    }
}


/* C4 exists iff two neighbours a,b of some vertex i have a common neighbour != i */
static int egf_has_C4_fast(void) {
    int i;
    for (i = 0; i < egf_n; i++) {
        unsigned long long bi = ~(1ULL << i);
        unsigned long long m0 = egf_nbrmask[egf_adj[i][0]];
        unsigned long long m1 = egf_nbrmask[egf_adj[i][1]];
        unsigned long long m2 = egf_nbrmask[egf_adj[i][2]];
        if ((m0 & m1 & bi) || (m0 & m2 & bi) || (m1 & m2 & bi)) return 1;
    }
    return 0;
}

/* 1 if a cycle of length exactly L exists */
static int egf_has_cycle(int L) {
    int r;
    if (L < 3 || L > egf_n) return 0;
    egf_L = L;
    for (r = 0; r + L <= egf_n; r++) {
        egf_root = r; egf_found = 0;
        if (L > 10) egf_bfs(r); else { int t; for (t = 0; t < egf_n; t++) egf_dist[t] = 0; }
        egf_used = 1ULL << r;
        egf_dfs(r, 0);
        if (egf_found) return 1;
    }
    return 0;
}

static void egf_to_graph6(char *out) {
    int n = egf_n, i, j, bitpos = 0, len, t;
    unsigned char bits[512];
    char *o = out;
    memset(bits, 0, sizeof(bits));
    for (j = 1; j < n; j++)
        for (i = 0; i < j; i++) {
            if ((egf_nbrmask[i] >> j) & 1ULL) bits[bitpos / 6] |= 1 << (5 - (bitpos % 6));
            bitpos++;
        }
    *o++ = (char)(n + 63);
    len = (bitpos + 5) / 6;
    for (t = 0; t < len; t++) *o++ = (char)(bits[t] + 63);
    *o = 0;
}


/* LEMMA UNDER TEST (crux of the truncation route):
     every cubic graph of girth 5 has a cycle of length 6, 7 or 8.
   We also measure the low-end spectrum so we can see HOW it fails if it does. */
static unsigned long long g5_seen=0, g5_girth5=0, g5_viol=0;
static unsigned long long g5_low[64];        /* histogram of least L>5 present */
static FILE *g5f=NULL;

int eg_filter(GRAPH g, int n){
    int i,j,L,least=0;
    g5_seen++;
    egf_n=n;
    for(i=0;i<n;i++){egf_nbrmask[i]=0;
        for(j=0;j<3;j++){egf_adj[i][j]=g[i][j];egf_nbrmask[i]|=1ULL<<g[i][j];}}
    if(!egf_has_cycle(5)) return 0;          /* need girth exactly 5 */
    g5_girth5++;
    for(L=6;L<=n && L<=12;L++) if(egf_has_cycle(L)){least=L;break;}
    if(least && least<64) g5_low[least]++;
    if(!(least>=6 && least<=8)){
        g5_viol++;
        if(!g5f){const char*p=getenv("EG_NEARFILE"); g5f=fopen(p?p:"g5viol.txt","a");}
        if(g5f){char b[512]; egf_to_graph6(b);
            fprintf(g5f,"%s n=%d leastAbove5=%d\n",b,n,least); fflush(g5f);}
        return 1;
    }
    return 0;
}

void eg_summary(void){
    int i;
    fprintf(stderr,"G5STAT seen=%llu girth5=%llu VIOLATIONS=%llu\n",g5_seen,g5_girth5,g5_viol);
    for(i=0;i<64;i++) if(g5_low[i]) fprintf(stderr,"G5STAT  least length >5 present = %d : %llu\n",i,g5_low[i]);
}
