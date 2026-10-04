"""vtx309.py -- verify: in a SUBCUBIC graph, a single-VERTEX odd cycle transversal forces
a single-EDGE one.  Hence tau(B) >= 2 implies no vertex meets all odd cycles.
Also measures, on bridges with tau >= 2, the vertex and edge transversal numbers and
nu_odd, to see what Lovasz-Slilaty would have to be applied to."""
import itertools, random, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split('"""', 2)[2])

def bipn(a):
    col={}
    for s in a:
        if s in col: continue
        col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in col: col[w]=1-col[u]; q.append(w)
                elif col[w]==col[u]: return None
    return col

def tau_vertex(E, cap=3):
    V=sorted({x for e in E for x in e}, key=str)
    for k in range(cap+1):
        for S in itertools.combinations(V,k):
            EE={e for e in E if not (set(e)&set(S))}
            if not EE: return k
            if bipn(adj_of(EE)) is not None: return k
    return cap+1

def tau_edge(E, cap=3):
    El=sorted(E,key=lambda s:sorted(map(str,s)))
    for k in range(cap+1):
        for S in itertools.combinations(El,k):
            if bipn(adj_of(set(El)-set(S))) is not None: return k
    return cap+1

# ---- exhaustive check of the lemma on small subcubic graphs -------------------
st=collections.Counter()
for n in range(4,8):
    allE=[frozenset(p) for p in itertools.combinations(range(n),2)]
    for k in range(n-1, min(len(allE), 3*n//2)+1):
        for S in itertools.combinations(allE,k):
            a=adj_of(set(S))
            if len(a)!=n: continue
            if max(len(a[v]) for v in a)>3: continue
            if bipn(a) is not None: continue         # bipartite: nothing to transverse
            tv=tau_vertex(set(S)); te=tau_edge(set(S))
            st['tested']+=1
            if tv==1:
                st['tau_vertex=1']+=1
                if te==1: st['  and tau_edge=1 (lemma holds)']+=1
                else: st['*** LEMMA FAILS ***']+=1
print('EXHAUSTIVE CHECK, subcubic non-bipartite graphs on 4..8 vertices')
for k in sorted(st): print('  %-36s %s'%(k,st[k]))

# ---- what tau >= 2 bridges look like ------------------------------------------
print()
print('BRIDGES WITH tau >= 2 (subdivided random cubic girth-5):')
rng=random.Random(17); pr=collections.Counter()
for _ in range(60):
    a,terms=make_bridge(rng.choice([10,12,14]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    te=tau_edge(E)
    if te<2: continue
    tv=tau_vertex(E)
    pr['tau_edge=%d, tau_vertex=%d'%(te,tv)]+=1
for k in sorted(pr): print('  %-30s %s'%(k,pr[k]))
