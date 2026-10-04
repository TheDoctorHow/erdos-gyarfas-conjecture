"""theta330.py -- (b) the transfer, and the theta reformulation.

B' = B + apex T* joined to b_1,b_2,b_3.  Claims to verify:
  1. B' is CUBIC - the terminals rise from degree 2 to 3 and T* has degree 3.
  2. In a subcubic graph every trail is a path, and edge-disjoint v-T* paths are
     automatically INTERNALLY VERTEX-DISJOINT (two through an interior vertex need 4 slots).
     So Ibrahimpur-Swamy transfers verbatim.
  3. Three internally disjoint ODD v-T* paths = a THETA subgraph with branch vertices v,T*
     and all three branches odd, equivalently all three of its cycles EVEN.
Measured: nu_odd(v) = max internally disjoint odd v-T* paths, over all v.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split(chr(34)*3,2)[2])

def all_simple(a, s, t, banned):
    out=[]; p=[s]; on=set(banned)|{s}
    def rec(u):
        for w in a[u]:
            if w==t: out.append(list(p)+[t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out

def nu_odd(a, v, t):
    """max internally disjoint ODD v-t paths, by brute force up to 3"""
    P=[p for p in all_simple(a,v,t,set()) if (len(p)-1)%2==1]
    best=0
    for r in (3,2,1):
        for combo in itertools.combinations(P, r):
            ok=True
            for x,y in itertools.combinations(combo,2):
                if set(x[1:-1]) & set(y[1:-1]): ok=False; break
            if ok: return r
    return best

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 233)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 10):
    a,terms=make_bridge(rng.choice([8,10]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    for x in terms: E.add(frozenset(('T*',x)))
    aB=adj_of(E)
    degs=collections.Counter(len(aB[u]) for u in aB)
    st['B-prime degree profile %s'%dict(degs)]+=1
    st['B-prime is CUBIC' if set(degs)=={3} else '*** B-prime NOT cubic ***']+=1
    best=0; nv=0
    for v in aB:
        if v=='T*': continue
        k=nu_odd(aB,v,'T*')
        best=max(best,k)
        if k>=3: nv+=1
    st['max nu_odd over v = %d'%best]+=1
    st['  vertices with nu_odd>=3: %d'%nv]+=1
for k in sorted(st): print('  %-42s %s'%(k,st[k]))
