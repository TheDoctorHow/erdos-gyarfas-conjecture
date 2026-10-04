/* count 0-mod-4 cycles, circumference, and pow2 mask for each graph on stdin */
#include "egcore.h"
int main(void){
    char line[8192]; EGGraph G;
    long long minc=-1; char worst[8192]={0};
    while(fgets(line,sizeof line,stdin)){
        char*p=strchr(line,'\n'); if(p)*p=0; if(!line[0])continue;
        if(eg_from_graph6(&G,line)){fprintf(stderr,"decode fail\n");return 1;}
        long long tot=0; int circ=0; unsigned mask=0; int L;
        for(L=3;L<=G.n;L++){
            long long r=eg_cycles_of_length(&G,L,0,2000000000LL);
            if(r<0){fprintf(stderr,"BUDGET %s L=%d\n",line,L);return 2;}
            if(r>0){circ=L; if((L&(L-1))==0) mask|=1u<<L;}
            if(r>0 && L%4==0) tot+=r;
        }
        printf("%s n=%d circ=%d n0mod4=%lld\n",line,G.n,circ,tot);
        if(minc<0||tot<minc){minc=tot; strcpy(worst,line);}
    }
    fprintf(stderr,"MIN #0mod4 cycles = %lld at %s\n",minc,worst);
    return 0;
}
