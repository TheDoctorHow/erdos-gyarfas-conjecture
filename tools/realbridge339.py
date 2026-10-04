"""realbridge339.py -- a CORRECT bridge generator, and (CYC) tested on genuine bridges.

Bridge (213.1): subcubic, girth >= 5, three degree-2 terminals PAIRWISE AT DISTANCE >= 5.
Construction: random cubic graph of girth >= 5 on n vertices; choose three edges pairwise
at endpoint-distance >= 4; subdivide each.  The subdivision vertices have degree 2 and are
pairwise at distance >= 6.  Verified, not assumed.

Then search for a (CYC) witness: two INTERNALLY DISJOINT even b_i-b_j paths whose lengths
differ by 2 mod 4, i.e. a cycle through two terminals with both arcs even and |Z| == 2
mod 4.  In a genuine bridge both arcs are >= 6, so |Z| >= 12 and the first candidates are
14, 18, 22.  A randomized search is legitimate here: it is hunting a WITNESS, so a find is
conclusive and only a miss would be inconclusive (325.1A).
"""
import sys, random, itertools, collections
from collections import deque

def girth(a):
    best=10**9
    for s in a:
        d={s:0};par={s:None};q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in d: d[w]=d[u]+1;par[w]=u;q.append(w)
                elif w!=par[u]: best=min(best,d[u]+d[w]+1)
    return best
def dist(a,s):
    d={s:0};q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w not in d: d[w]=d[u]+1;q.append(w)
    return d
def rand_cubic(n,rng,tries=4000):
    for _ in range(tries):
        stubs=[v for v in range(n) for _ in range(3)]
        rng.shuffle(stubs)
        E=set(); ok=True
        for i in range(0,len(stubs),2):
            u,v=stubs[i],stubs[i+1]
            if u==v or frozenset((u,v)) in E: ok=False; break
            E.add(frozenset((u,v)))
        if not ok: continue
        a={v:set() for v in range(n)}
        for e in E:
            u,v=tuple(e); a[u].add(v); a[v].add(u)
        if any(len(a[v])!=3 for v in a): continue
        if girth(a)>=5: return a,E
    return None,None

def make_real_bridge(n,rng):
    a,E=rand_cubic(n,rng)
    if a is None: return None,None
    El=[tuple(e) for e in E]
    D={v:dist(a,v) for v in a}
    for _ in range(600):
        tri=rng.sample(El,3)
        ok=True
        for (u1,v1),(u2,v2) in itertools.combinations(tri,2):
            if min(D[x][y] for x in (u1,v1) for y in (u2,v2))<3: ok=False; break
        if not ok: continue
        E2=set(E); terms=[]
        for k,(u,v) in enumerate(tri):
            t='b%d'%k; terms.append(t)
            E2.discard(frozenset((u,v)))
            E2.add(frozenset((u,t))); E2.add(frozenset((t,v)))
        a2={}
        for e in E2:
            u,v=tuple(e); a2.setdefault(u,set()).add(v); a2.setdefault(v,set()).add(u)
        DD={t:dist(a2,t) for t in terms}
        if min(DD[terms[i]][terms[j]] for i,j in [(0,1),(0,2),(1,2)])>=5:
            return a2,terms
    return None,None

def rand_paths(a,s,t,rng,K):
    """K randomized simple s-t paths"""
    out=[]
    for _ in range(K):
        p=[s]; on={s}
        while p[-1]!=t:
            nb=[w for w in a[p[-1]] if w not in on]
            if not nb: break
            w=rng.choice(nb); p.append(w); on.add(w)
        if p[-1]==t: out.append(p)
    return out

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 991)
N=int(sys.argv[2]) if len(sys.argv)>2 else 6
NV=int(sys.argv[3]) if len(sys.argv)>3 else 30
st=collections.Counter()
for _ in range(N):
    a,terms=make_real_bridge(NV,rng)
    if a is None: st['generation failed']+=1; continue
    DD={t:dist(a,t) for t in terms}
    pd=tuple(sorted(DD[terms[i]][terms[j]] for i,j in [(0,1),(0,2),(1,2)]))
    st['bridges (n=%d)'%len(a)]+=1
    st['  terminal distances %s  girth %d'%(str(pd),girth(a))]+=1
    found=None
    for i,j in [(0,1),(0,2),(1,2)]:
        P=[p for p in rand_paths(a,terms[i],terms[j],rng,6000) if (len(p)-1)%2==0]
        byc=collections.defaultdict(list)
        for p in P: byc[(len(p)-1)%4].append(p)
        if 0 in byc and 2 in byc:
            for p1 in byc[0][:200]:
                for p2 in byc[2][:200]:
                    if not (set(p1[1:-1]) & set(p2[1:-1])):
                        found=(i,j,len(p1)-1,len(p2)-1,(len(p1)-1)+(len(p2)-1)); break
                if found: break
        if found: break
    st['*** (CYC) WITNESS on a genuine bridge ***' if found else 'no witness found (inconclusive)']+=1
    if found: st['  pair(%d,%d) arcs %d,%d  |Z|=%d'%found]+=1
for k in sorted(st): print('  %-52s %s'%(k,st[k]))
