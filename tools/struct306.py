"""struct306.py -- what does 'the collapse fails' force on H?

If for EVERY chordless x-b_2 path P every component of R = H - V(P) carries at most TWO
landing points, then each such component is joined to P by exactly two edges: a 2-EDGE-CUT
of H.  This inspects the component structure of R in the 84 resistant instances.
"""
import itertools, collections
from collections import deque
exec(open('tools/both306.py').read().split('st=collections.Counter()')[0].split('"""',2)[2])

st=collections.Counter(); prof=collections.Counter()
E0=heawood(); El=sorted(E0,key=lambda f:sorted(map(str,f)))
for e1,e2 in itertools.combinations(El,2):
    if e1&e2: continue
    H=set(E0)-{e1,e2}; aH=adj(H)
    if not two_conn(aH) or girth(aH)<6: continue
    cH=bip(aH)
    if cH is None: continue
    d2=[v for v in aH if len(aH[v])==2]
    if len(d2)!=4: continue
    X=[v for v in d2 if cH[v]==0]; Y=[v for v in d2 if cH[v]==1]
    if len(X)!=2 or len(Y)!=2: continue
    for x in X:
        b3=[v for v in X if v!=x][0]
        for b2 in Y:
            y=[v for v in Y if v!=b2][0]
            if collapse_fires(aH,x,b2): continue
            st['resistant instances']+=1
            # profile the component structure over all chordless x-b2 paths
            maxland=0; ncut2=0; npaths=0
            for p in paths(aH,x,b2):
                Pv=set(p); ok=True
                for i,v in enumerate(p):
                    inn=[w for w in aH[v] if w in Pv]
                    if len(inn)>(1 if i in (0,len(p)-1) else 2): ok=False;break
                if not ok: continue
                npaths+=1
                land={}
                for i,v in enumerate(p):
                    off=[w for w in aH[v] if w not in Pv]
                    if len(off)==1: land[i]=off[0]
                R=set(aH)-Pv
                if not R: continue
                for C in comps(sub(aH,R)):
                    k=len([i for i in land if land[i] in C])
                    maxland=max(maxland,k)
                    # edges from C to P
                    cut=sum(1 for v in C for w in aH[v] if w in Pv)
                    if cut==2: ncut2+=1
                    prof['component: %d landing points, %d edges to P'%(k,cut)]+=1
            st['max landing points in any component = %d'%maxland]+=1
            st['chordless paths examined']+=npaths
for k in sorted(st): print('  %-46s %s'%(k,st[k]))
print()
for k in sorted(prof): print('   %-46s %s'%(k,prof[k]))
