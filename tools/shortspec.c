/* shortspec.c -- print which cycle lengths 3..L are present (exact). */
#include "egcore.h"
int main(int argc,char**argv){
    char line[1<<17]; int L=(argc>1)?atoi(argv[1]):12;
    while(fgets(line,sizeof(line),stdin)){
        EGGraph G; int k,mind; char*p=strchr(line,'\n'); if(p)*p=0; p=strchr(line,'\r'); if(p)*p=0;
        if(!line[0])continue;
        if(eg_from_graph6(&G,line)!=0){fprintf(stderr,"BAD\n");continue;}
        mind=99; for(k=0;k<G.n;k++) if(G.deg[k]<mind) mind=G.deg[k];
        printf("n=%d mindeg=%d  present:",G.n,mind);
        for(k=3;k<=L&&k<=G.n;k++)
            if(eg_cycles_of_length(&G,k,1,2000000000LL)>0) printf(" %d",k);
        printf("\n");
    }
    return 0;
}
