"""faneven.py -- structure behind (FAN-EVEN).

For each (N3) graph H and pair (x,b_2): over all 6-cycles Z, record
  * the EXITS of Z (vertices of Z with an edge off Z),
  * how many exits, and how many SAME-CLASS pairs of exits,
  * which arc parities are actually achieved by x-b_2 paths meeting Z in a single arc.
Question: is "Z has two same-class exits" enough to force an even arc?
"""
import itertools, collections
from collections import deque
exec(open('tools/n3search.py').read().split('def heawood()')[0].split('"""',2)[2])
def heawood():
    E=set()
    for i in range(14): E.add(frozenset((i,(i+1)%14)))
    for i in range(0,14,2): E.add(frozenset((i,(i+5)%14)))
    return E
def eset(p): return frozenset(frozenset((p[i],p[i+1])) for i in range(len(p)-1))

def six_cycles(a):
    out=set(); vs=sorted(a,key=str); idx={v:i for i,v in enumerate(vs)}
    for i,s in enumerate(vs):
        p=[s]; on={s}
        def rec(u):
            if len(p)>6: return
            for w in a[u]:
                if idx[w]<i: continue
                if w==s and len(p)==6: out.add(frozenset(p))
                elif w not in on and len(p)<6:
                    on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
        rec(s)
    return out

st=collections.Counter(); corr=collections.Counter()
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
    SC=six_cycles(aH)
    for x in X:
        for b2 in Y:
            P=[p for p in paths(aH,x,b2)]
            st['pairs']+=1
            good=0
            for Zv in SC:
                Zed={frozenset((u,w)) for u in Zv for w in aH[u] if w in Zv}
                if len(Zed)!=6: continue
                exits=[v for v in Zv if any(w not in Zv for w in aH[v])]
                same=[(u,v) for u,v in itertools.combinations(exits,2) if cH[u]==cH[v]]
                pars=set()
                for p in P:
                    common=[v for v in p if v in Zv]
                    ed=eset(p)&Zed
                    if not ed: continue
                    # single arc test: the common vertices form a contiguous run on p and
                    # induce a path in Z with |ed| edges and |ed|+1 vertices
                    if len(common)!=len(ed)+1: continue
                    ai=[i for i,v in enumerate(p) if v in Zv]
                    if ai[-1]-ai[0]!=len(ai)-1: continue
                    pars.add(len(ed)%2)
                if 0 in pars: good+=1
                corr['exits=%d same-class-exit-pairs=%d -> even arc achievable: %s'
                     % (len(exits),len(same),0 in pars)]+=1
            st['pairs with an even-arc 6-cycle' if good else
               '*** NO even-arc 6-cycle ***']+=1
for k in sorted(st): print('  %-40s %s' % (k,st[k]))
print()
for k in sorted(corr): print('   %-62s %s' % (k,corr[k]))
