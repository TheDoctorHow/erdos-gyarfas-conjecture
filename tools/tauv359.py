"""tauv359.py -- tau_v, the minimum odd-cycle VERTEX transversal, in genuine bridges.

(NOT-ALL-TRANSVERSAL) asks that some shortest even terminal cycle not be an odd-cycle
transversal.  A hexagon has SIX vertices, so:
        tau_v >= 7   ==>   NO hexagon can be a transversal   ==>   gap (B) holds for hexagons.
And 256.1 caps tau_v <= 7 in a surviving bridge.  So the whole question is whether tau_v sits
at 7 or drops to 6 or below.  Measured directly: smallest vertex set whose deletion leaves a
bipartite graph.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def bip_after(a, S):
    col={}
    for s in a:
        if s in S or s in col: continue
        col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w in S: continue
                if w not in col: col[w]=1-col[u]; q.append(w)
                elif col[w]==col[u]: return False
    return True

def tau_v(a, cap=8):
    V=sorted(a, key=str)
    for k in range(cap+1):
        for S in itertools.combinations(V,k):
            if bip_after(a,set(S)): return k,S
    return cap+1,None

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 6101)
st=collections.Counter()
for NV in (30,40):
    for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 3):
        a,terms=make_real_bridge(NV,rng)
        if a is None: continue
        k,S=tau_v(a,cap=6)          # only need to know whether it is <= 6
        if S is None:
            st['n=%-3d  tau_v >= 7   -> NO hexagon can be a transversal'%len(a)]+=1
        else:
            st['n=%-3d  tau_v = %d   (<= 6: a hexagon COULD be one)'%(len(a),k)]+=1
for k in sorted(st): print('  %s  %s'%(k.ljust(56),st[k]))
