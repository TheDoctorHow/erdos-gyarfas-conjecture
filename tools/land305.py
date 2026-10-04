"""land305.py -- is the LANDING-POINT CONGRUENCE itself always violated?

Rigidity forces, for any chordless x-b_2 path P with landing points p_i, and any i<j:
        every p_i-p_j path of R = H - V(P)  has length == j - i - 2  (mod 4).
Rigidity is false in all 672 instances, so some (P,i,j,path) must violate it.  This finds
the violation and reports how it arises: by a TRIPLE CONCATENATION (588 cases, the elegant
mechanism) or otherwise.
"""
import itertools, collections
from collections import deque
exec(open('tools/n3search.py').read().split('def heawood()')[0].split('"""',2)[2])
def heawood():
    E=set()
    for i in range(14): E.add(frozenset((i,(i+1)%14)))
    for i in range(0,14,2): E.add(frozenset((i,(i+5)%14)))
    return E
def sub(a,S): return {v:{w for w in a[v] if w in S} for v in S}

st=collections.Counter(); how=collections.Counter()
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
        for b2 in Y:
            st['pairs']+=1
            viol=None
            for p in paths(aH,x,b2):
                Pv=set(p); ok=True
                for i,v in enumerate(p):
                    inn=[w for w in aH[v] if w in Pv]
                    if len(inn)>(1 if i in (0,len(p)-1) else 2): ok=False;break
                if not ok: continue
                land={}
                for i,v in enumerate(p):
                    off=[w for w in aH[v] if w not in Pv]
                    if len(off)==1: land[i]=off[0]
                R=set(aH)-Pv
                if len(R)<2: continue
                aR=sub(aH,R)
                idxs=sorted(land)
                for i,j in itertools.combinations(idxs,2):
                    pi,pj=land[i],land[j]
                    if pi==pj: continue
                    for q in paths(aR,pi,pj):
                        if (len(q)-1-(j-i-2))%4 != 0:
                            viol=(len(p)-1,i,j,len(q)-1); break
                    if viol: break
                if viol: break
            if viol:
                st['landing congruence VIOLATED']+=1
                how['|P|=%d gap=%d detour=%d'%(viol[0],viol[2]-viol[1],viol[3])]+=1
            else:
                st['*** no violation found ***']+=1
for k in sorted(st): print('  %-42s %s'%(k,st[k]))
print('  (sample of how the first violation looks)')
for k in sorted(how)[:8]: print('    %-28s %s'%(k,how[k]))
