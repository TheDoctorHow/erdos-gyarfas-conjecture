/* gadget.c -- enumerate 1-terminal gadgets (one vertex of degree 2, rest degree 3)
 * from cubic base graphs, and test them for power-of-two cycles.
 *
 * WHY.  If a cubic graph G has a bridge, no cycle uses it, so the cycle spectrum
 * of G is the union of the spectra of the two sides.  A LEAF block of the bridge
 * tree of a cubic graph is exactly a "1-terminal gadget" H: one vertex of degree
 * 2, all others degree 3.  Conversely, joining the terminals of two gadgets by an
 * edge gives a cubic graph whose spectrum is the union of the two spectra.
 * Hence:  a cubic counterexample with a bridge exists  <=>  a 1-terminal gadget
 * with no cycle of length a power of two exists.
 * A gadget of order k has circumference <= k, so for k <= 31 it only has to dodge
 * {4,8,16}: every larger power of two is STRUCTURALLY absent.  Doubling such a
 * gadget yields a cubic counterexample on 2k <= 62 vertices.
 *
 * COMPLETENESS.  Let t be the terminal, with neighbours p,q.
 *   - p not adjacent to q:  M := H - t + pq is cubic on k-1 vertices and
 *     H = subdivide(M, pq).                                   [FAMILY 1]
 *   - p adjacent to q:  t,p,q is a triangle; let x,y be the third neighbours of
 *     p,q.  M := H - {t,p,q} + xy is cubic on k-3 vertices.    [FAMILY 2]
 * So every 1-terminal gadget arises from a cubic graph by one of the two
 * operations below.
 *
 * Usage:  snarkhunter <m> 3 o g s | gadget <1|2>
 */
#include "egcore.h"

static long long n_base=0,n_g=0,n_noC4=0,n_noC4C8=0,n_ctr=0;
static long long histo[8];
static FILE *nearf=NULL;

static void test(EGGraph *H,int k){
    unsigned mask,unk; int j,low=-1; char g6[1024];
    n_g++;
    if(eg_has_C4(H))return;
    n_noC4++;
    if(eg_cycles_of_length(H,8,1,200000000LL)>0)return;
    n_noC4C8++;
    mask=eg_pow2_mask(H,800000000LL,&unk);
    for(j=2;j<32;j++) if(mask&(1u<<j)){low=j;break;}
    if(low>=0&&low<8)histo[low]++;
    eg_to_graph6(H,g6);
    if(nearf){fprintf(nearf,"%s k=%d mask=%u unk=%u lowest=%d\n",g6,k,mask,unk,low>=0?(1<<low):-1);fflush(nearf);}
    if(mask==0&&unk==0){
        n_ctr++;
        printf("GADGET_COUNTEREXAMPLE %s k=%d  -> double it for a cubic counterexample on %d vertices\n",g6,k,2*k);
        fflush(stdout);
        fprintf(stderr,"\n*** GADGET COUNTEREXAMPLE *** %s k=%d\n",g6,k);
    }
}

int main(int argc,char**argv){
    char line[1<<16]; int fam=(argc>1)?atoi(argv[1]):1; int i;
    const char*nf=getenv("EG_NEARFILE"); if(nf)nearf=fopen(nf,"a");
    while(fgets(line,sizeof(line),stdin)){
        EGGraph M,H; char*p=strchr(line,'\n'); if(p)*p=0; p=strchr(line,'\r'); if(p)*p=0;
        if(!line[0])continue;
        if(eg_from_graph6(&M,line)!=0){fprintf(stderr,"BAD %s\n",line);continue;}
        n_base++;
        {
            int m=M.n,u,a,b,c,d;
            /* PRUNE: subdividing/replacing edge uv removes only the C4s of M that
               USE uv.  So the gadget is C4-free only if every C4 of M contains uv.
               Collect the intersection of the edge sets of all C4s of M. */
            unsigned long long okmask = ~0ULL;      /* bit = index into eid[][] */
            int eid[EG_MAXV][EG_MAXV], ne=0;
            for(c=0;c<m;c++)for(d=0;d<m;d++) eid[c][d]=-1;
            for(c=0;c<m;c++)for(d=0;d<M.deg[c];d++){int w=M.adj[c][d]; if(w<c)continue;
                eid[c][w]=eid[w][c]=ne++;}
            for(c=0;c<m && okmask;c++)
              for(a=0;a<M.deg[c];a++)
                for(b=a+1;b<M.deg[c];b++){
                  int p1=M.adj[c][a],p2=M.adj[c][b],y;
                  for(y=0;y<M.deg[p1];y++){int w=M.adj[p1][y];
                    if(w==c||!M.is_adj[w][p2])continue; if(w<c)continue;
                    okmask &= (1ULL<<eid[c][p1])|(1ULL<<eid[p1][w])|
                              (1ULL<<eid[w][p2])|(1ULL<<eid[p2][c]);}
                }
            if(!okmask) continue;
            for(u=0;u<m;u++)for(a=0;a<M.deg[u];a++){
                int v=M.adj[u][a]; if(v<u)continue;      /* each edge uv once */
                if(!((okmask>>eid[u][v])&1ULL)) continue;
                if(fam==1){
                    /* subdivide uv: new vertex m is the terminal */
                    eg_clear(&H,m+1);
                    for(c=0;c<m;c++)for(d=0;d<M.deg[c];d++){int w=M.adj[c][d];
                        if(w<c)continue; if((c==u&&w==v)||(c==v&&w==u))continue;
                        eg_add_edge(&H,c,w);}
                    eg_add_edge(&H,u,m); eg_add_edge(&H,m,v);
                    test(&H,m+1);
                }else{
                    /* replace uv by u-A-B-v and add terminal T joined to A,B */
                    int A=m,B=m+1,T=m+2;
                    eg_clear(&H,m+3);
                    for(c=0;c<m;c++)for(d=0;d<M.deg[c];d++){int w=M.adj[c][d];
                        if(w<c)continue; if((c==u&&w==v)||(c==v&&w==u))continue;
                        eg_add_edge(&H,c,w);}
                    eg_add_edge(&H,u,A); eg_add_edge(&H,A,B); eg_add_edge(&H,B,v);
                    eg_add_edge(&H,T,A); eg_add_edge(&H,T,B);
                    test(&H,m+3);
                }
            }
        }
    }
    fprintf(stderr,"GAD fam=%d bases=%lld gadgets=%lld noC4=%lld noC4C8=%lld ctr=%lld\n",
            fam,n_base,n_g,n_noC4,n_noC4C8,n_ctr);
    for(i=0;i<8;i++) if(histo[i]) fprintf(stderr,"GAD   noC4C8 smallest 2^k=%d : %lld\n",1<<i,histo[i]);
    if(nearf)fclose(nearf);
    return 0;
}
