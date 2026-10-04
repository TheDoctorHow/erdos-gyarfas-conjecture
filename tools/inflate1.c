/* inflate1.c -- COMPLETE triangle reduction for the Erdos-Gyarfas hunt.
 *
 * THEOREM.  Let G be a connected cubic C4-free graph containing a triangle
 * T = {x,y,z}.  Then the three outside neighbours of x,y,z are distinct (else
 * G has a C4), so contracting T to a single vertex v yields a connected cubic
 * simple graph H on n-2 vertices, and G = inflate(H, v).  Moreover
 *
 *    inflate(H,v) is C4-free   <=>   (a) no triangle of H passes through v, and
 *                                    (b) every C4 of H passes through v.
 *
 * Proof of <=>.  A cycle of G meets the triangle T in 0 or >= 2 vertices (a
 * cycle through x_i uses two of its three edges, at least one inside T).
 *  - C4s of G avoiding T are exactly C4s of H avoiding v          -> (b)
 *  - a C4 meeting T in 2 vertices is x_i x_j b a x_i where a,b are the (unique,
 *    distinct) outside neighbours of x_i,x_j and a~b; in H that is the triangle
 *    v a b                                                        -> (a)
 *  - a C4 meeting T in 3 vertices would need one outside vertex adjacent to two
 *    triangle vertices: impossible.
 * Hence every connected cubic C4-free graph on n vertices either has girth >= 5,
 * or is inflate(H,v) for some cubic H on n-2 vertices and some vertex v
 * satisfying (a),(b).  This is COMPLETE (my first attempt, which contracted all
 * triangles at once and assumed the base was triangle-free, was NOT: a base
 * triangle whose vertices are 2 or 3 contracted vertices creates no C4).
 *
 * Usage:  snarkhunter <n-2> 3 o g s | inflate1
 */
#include "egcore.h"

static long long n_base=0, n_cand=0, n_noC4=0, n_noC4C8=0, n_ctr=0, n_skiptri=0;
static long long histo[8];
static FILE *nearf=NULL;

/* build inflate(H,v) */
static void build(const EGGraph *H, int v, EGGraph *G) {
    int m=H->n, i, p, idx=0, img[EG_MAXV], base;
    for (i=0;i<m;i++) img[i] = (i==v) ? -1 : idx++;
    base = idx; idx += 3;
    eg_clear(G, idx);
    eg_add_edge(G, base+0, base+1);
    eg_add_edge(G, base+1, base+2);
    eg_add_edge(G, base+2, base+0);
    for (i=0;i<m;i++) for (p=0;p<H->deg[i];p++) {
        int j=H->adj[i][p], q, a, b;
        if (j<i) continue;
        for (q=0;q<H->deg[j];q++) if (H->adj[j][q]==i) break;
        a = (i==v) ? base+p : img[i];
        b = (j==v) ? base+q : img[j];
        eg_add_edge(G, a, b);
    }
}

int main(int argc, char **argv) {
    char line[1<<16];
    const char *nf=getenv("EG_NEARFILE");
    int i;
    if (nf) nearf=fopen(nf,"a");
    while (fgets(line,sizeof(line),stdin)) {
        EGGraph H, G;
        char *p=strchr(line,'\n'); if(p)*p=0; p=strchr(line,'\r'); if(p)*p=0;
        if(!line[0]) continue;
        if (eg_from_graph6(&H,line)!=0){fprintf(stderr,"BAD %s\n",line);continue;}
        n_base++;
        {
            int m=H.n, v, u, a, b, x;
            /* candidate v must lie on EVERY C4 of H; collect intersection */
            unsigned long long cand = (m>=64)?~0ULL:((1ULL<<m)-1);
            int anyc4=0;
            for (u=0;u<m && cand;u++)
              for (a=0;a<H.deg[u];a++)
                for (b=a+1;b<H.deg[u];b++){
                  int p1=H.adj[u][a], p2=H.adj[u][b], y;
                  for (y=0;y<H.deg[p1];y++){
                    int w=H.adj[p1][y];
                    if(w==u||!H.is_adj[w][p2])continue;
                    if(w<u)continue;
                    anyc4=1;
                    cand &= (1ULL<<u)|(1ULL<<p1)|(1ULL<<w)|(1ULL<<p2);
                  }
                }
            (void)anyc4;
            for (v=0; v<m; v++) {
                if(!((cand>>v)&1ULL)) continue;
                /* (a) no triangle through v */
                { int t=0,q,r;
                  for(q=0;q<H.deg[v];q++) for(r=q+1;r<H.deg[v];r++)
                      if(H.is_adj[H.adj[v][q]][H.adj[v][r]]) t=1;
                  if(t){n_skiptri++;continue;} }
                n_cand++;
                build(&H,v,&G);
                if (eg_has_C4(&G)) { fprintf(stderr,"INTERNAL: C4 leaked\n"); continue; }
                n_noC4++;
                if (eg_cycles_of_length(&G,8,1,100000000LL)>0) continue;
                n_noC4C8++;
                { unsigned mask,unk; int k,low=-1; char g6[512];
                  mask=eg_pow2_mask(&G,400000000LL,&unk);
                  for(k=2;k<32;k++) if(mask&(1u<<k)){low=k;break;}
                  if(low>=0&&low<8) histo[low]++;
                  eg_to_graph6(&G,g6);
                  if(nearf){fprintf(nearf,"%s n=%d mask=%u unk=%u lowest=%d\n",g6,G.n,mask,unk,low>=0?(1<<low):-1);fflush(nearf);}
                  if(mask==0&&unk==0){n_ctr++;printf("COUNTEREXAMPLE %s n=%d\n",g6,G.n);fflush(stdout);
                    fprintf(stderr,"\n*** COUNTEREXAMPLE *** %s n=%d\n",g6,G.n);}
                }
            }
        }
    }
    fprintf(stderr,"INF1 bases=%lld cand=%lld skiptri=%lld noC4=%lld noC4C8=%lld ctr=%lld\n",
            n_base,n_cand,n_skiptri,n_noC4,n_noC4C8,n_ctr);
    for(i=0;i<8;i++) if(histo[i]) fprintf(stderr,"INF1   noC4C8 smallest 2^k=%d : %lld\n",1<<i,histo[i]);
    if(nearf)fclose(nearf);
    (void)argc;(void)argv;
    return 0;
}
