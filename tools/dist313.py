"""dist313.py -- can a RIGID STALL BLOCK occur at tau = 2, using only PROVED facts?

At tau = 2 the ambient degree-2 vertices of W = B - N are the 3 terminals plus the 4
endpoints of the two N-edges.  Proved (310.2):
   * terminals are pairwise at distance >= 5,
   * the two endpoints of ONE N-edge are at distance >= 4.
A block K of the descent has s(K) degree-2 vertices, two of which are its ends a,b; the
other s-2 must be ambient degree-2 vertices lying inside K, and distances inside K are at
least the ambient ones.  So a rigid stall block with s <= 7 can occur only if its s-2
non-end degree-2 vertices can be labelled as (at most one terminal) plus N-endpoints
forming pairs at mutual distance >= 4.
This tests every rigid ladder found by 312 with s <= 7.
"""
import itertools, collections
from collections import deque
exec(open('tools/ladder312.py').read().split('st=collections.Counter()')[0].split('"""',2)[2])

def dists(a,s):
    d={s:0};q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w not in d: d[w]=d[u]+1;q.append(w)
    return d

st=collections.Counter(); shown=0
for L in range(4, 11):
    V=list(range(L+1))
    for M in matchings(V):
        if not M: continue
        ranges=[];ok=True
        for (p,q) in M:
            sp=q-p
            opts=[l for l in range(1,5) if (l-sp)%2==0 and l+sp>=6]
            if not opts: ok=False;break
            ranges.append(opts)
        if not ok: continue
        for combo in itertools.product(*ranges):
            E=set()
            for i in range(L): E.add(frozenset(('v%d'%i,'v%d'%(i+1))))
            idc=0
            for (p,q),l in zip(M,combo):
                prev='v%d'%p
                for k in range(l-1):
                    nv='c%d_%d'%(idc,k); E.add(frozenset((prev,nv))); prev=nv
                E.add(frozenset((prev,'v%d'%q))); idc+=1
            a=adj(E)
            if max(len(a[v]) for v in a)>3: continue
            if girth(a)<6 or not two_conn(a): continue
            Ls={l for l in paths(a,'v0','v%d'%L)}
            if len({x%4 for x in Ls})!=1: continue
            d2=[v for v in a if len(a[v])==2]
            s=len(d2)
            if s>7: continue
            st['rigid stall blocks with s<=7']+=1
            ends={'v0','v%d'%L}
            non=[v for v in d2 if v not in ends]
            # try: choose <=1 of them to be a terminal, rest are N-endpoints in pairs >=4
            DD={v:dists(a,v) for v in non}
            feasible=False
            for term in [None]+non:
                rest=[v for v in non if v!=term]
                if len(rest)%2: continue
                if len(rest)>4: continue          # tau=2 gives only 4 N-endpoints
                # perfect matching of rest with all pairs at distance >= 4
                def go(rem):
                    if not rem: return True
                    x=rem[0]
                    for y in rem[1:]:
                        if DD[x].get(y,99)>=4 and go([z for z in rem[1:] if z!=y]): return True
                    return False
                if go(rest): feasible=True;break
            st['  CAN occur at tau=2' if feasible else '  EXCLUDED at tau=2']+=1
            if not feasible and shown<5:
                shown+=1
                print('   excluded: L=%d chords=%s len=%s  s=%d n=%d'%(L,M,combo,s,len(a)))
for k in sorted(st): print('  %-34s %s'%(k,st[k]))
