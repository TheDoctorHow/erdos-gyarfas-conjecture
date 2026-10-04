"""caveat333.py -- does the 2-connectivity caveat of 330.1 ever bite?

330.3 proves tau(B') >= 6 forces an all-even fan, ASSUMING B' - F_v stays 2-connected for
the 5-edge certificate F_v.  B' is cubic and 3-connected, hence 3-EDGE-connected.  Two
questions:
  (1) how many components can 5 edge-deletions create?  Each component of B'-F has edge
      boundary >= 3, and each deleted edge serves at most 2 boundaries, so 3c <= 2|F| = 10,
      giving c <= 3.  Verified empirically.
  (2) how often does deleting 5 edges destroy 2-CONNECTIVITY (create a cut vertex or
      disconnect)?  If it is rare the caveat is cosmetic; if common it is essential.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split(chr(34)*3,2)[2])

def comps(a, vs):
    seen=set(); c=0
    for s in vs:
        if s in seen: continue
        c+=1; seen.add(s); q=deque([s])
        while q:
            u=q.popleft()
            for w in a.get(u,()):
                if w not in seen: seen.add(w); q.append(w)
    return c
def two_conn(a, vs):
    if comps(a,vs)>1: return False
    for x in vs:
        rest=[v for v in vs if v!=x]
        if not rest: continue
        seen={rest[0]}; q=deque([rest[0]])
        while q:
            u=q.popleft()
            for w in a.get(u,()):
                if w!=x and w not in seen: seen.add(w); q.append(w)
        if len(seen)!=len(rest): return False
    return True

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 271)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 6):
    a,terms=make_bridge(10,rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    for x in terms: E.add(frozenset(('T*',x)))
    V=sorted(set().union(*[set(e) for e in E]), key=str)
    El=sorted(E,key=lambda s:sorted(map(str,s)))
    st['bridges']+=1
    st['  |V(B-prime)| = %d, |E| = %d'%(len(V),len(El))]+=1
    tot=0; broke=0; mx=1
    for F in itertools.combinations(El,5):
        aF=adj_of(set(El)-set(F))
        tot+=1
        c=comps(aF,V); mx=max(mx,c)
        if not two_conn(aF,V): broke+=1
    st['  5-edge deletions tested: %d'%tot]+=1
    st['  max components created: %d'%mx]+=1
    st['  fraction destroying 2-connectivity: %.2f'%(broke/max(tot,1))]+=1
for k in sorted(st): print('  %-46s %s'%(k,st[k]))
