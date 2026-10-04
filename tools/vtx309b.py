"""vtx309b.py -- fast exhaustive check of the Vertex-to-Edge Transversal Lemma:
in a SUBCUBIC graph, if some vertex meets every odd cycle then some edge does."""
import itertools, collections
from collections import deque
def adj(E):
    a={}
    for e in E:
        u,v=tuple(e); a.setdefault(u,set()).add(v); a.setdefault(v,set()).add(u)
    return a
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
st=collections.Counter()
for n in range(4,9):
    allE=[frozenset(p) for p in itertools.combinations(range(n),2)]
    for k in range(n-1, 3*n//2+1):
        for S in itertools.combinations(allE,k):
            E=set(S); a=adj(E)
            if len(a)!=n: continue
            if max(len(a[v]) for v in a)>3: continue
            if bipn(a) is not None: continue
            # tau_vertex == 1 ?
            tv1=[v for v in a if bipn(adj({e for e in E if v not in e})) is not None]
            if not tv1: continue
            st['graphs with a single-vertex transversal']+=1
            te1=[e for e in E if bipn(adj(E-{e})) is not None]
            if te1: st['  ... also a single-edge transversal']+=1
            else:
                st['*** LEMMA FAILS ***']+=1
                print('FAIL n=%d E=%s v=%s'%(n,sorted(tuple(sorted(e)) for e in E),tv1))
            # sharper: for a degree-3 transversal vertex, is the edge unique+canonical?
            for v in tv1:
                if len(a[v])==3:
                    inc=[e for e in te1 if v in e]
                    st['deg-3 transversal vertex, incident bipartizing edges: %d'%len(inc)]+=1
for k in sorted(st): print('  %-52s %s'%(k,st[k]))
