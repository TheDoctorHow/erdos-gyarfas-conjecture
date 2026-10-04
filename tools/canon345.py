"""canon345.py -- what distinguishes the clean meetings where the arc-shift fires?

If a CANONICAL choice - smallest even cycle D, shortest tail Q_2 - always satisfies the
construction's hypothesis, then proving (ODD-DISJ') is a matter of making that choice.
Correlate success against |D| and |Q_2|.  Self-contained.
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
    El=[tuple(e) for e in E]; D={v:dist(a,v) for v in a}
    for _ in range(600):
        tri=rng.sample(El,3); ok=True
        for (u1,v1),(u2,v2) in itertools.combinations(tri,2):
            if min(D[x][y] for x in (u1,v1) for y in (u2,v2))<3: ok=False; break
        if not ok: continue
        E2=set(E); terms=[]
        for k,(u,v) in enumerate(tri):
            t='b%d'%k; terms.append(t)
            E2.discard(frozenset((u,v))); E2.add(frozenset((u,t))); E2.add(frozenset((t,v)))
        a2={}
        for e in E2:
            u,v=tuple(e); a2.setdefault(u,set()).add(v); a2.setdefault(v,set()).add(u)
        DD={t:dist(a2,t) for t in terms}
        if min(DD[terms[i]][terms[j]] for i,j in [(0,1),(0,2),(1,2)])>=5: return a2,terms
    return None,None
def rand_paths(a,s,t,rng,K):
    out=[]
    for _ in range(K):
        p=[s]; on={s}
        while p[-1]!=t:
            nb=[w for w in a[p[-1]] if w not in on]
            if not nb: break
            w=rng.choice(nb); p.append(w); on.add(w)
        if p[-1]==t: out.append(p)
    return out
def rand_cycle(a,rng):
    """returns the actual CYCLE as a vertex list; len == cycle length"""
    u=rng.choice(list(a)); v=rng.choice(list(a[u])); p=[u]; on={u,v}
    while True:
        nb=[w for w in a[p[-1]] if w not in on]
        if not nb: return None
        w=rng.choice(nb); p.append(w); on.add(w)
        if v in a[w] and len(p)>=3: return p+[v]
        if len(p)>len(a): return None
def arc_span(D,S):
    n=len(D); idx=[i for i,x in enumerate(D) if x in S]
    if len(idx)!=len(S) or not idx: return None
    for s in idx:
        run=[D[(s+k)%n] for k in range(len(S))]
        if set(run)==S: return run
    return None
def comp_of(a,drop,s):
    if s in drop: return None
    seen={s}; q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w in drop or w in seen: continue
            seen.add(w); q.append(w)
    return seen
def nonbip(a,vs):
    col={}
    for s in vs:
        if s in col: continue
        col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in vs: continue
                if w not in col: col[w]=1-col[u]; q.append(w)
                elif col[w]==col[u]: return True
    return False
def two_conn(a,vs):
    vs=set(vs)
    if len(vs)<3: return False
    for x in vs:
        rest=vs-{x}
        c=comp_of(a,(set(a)-vs)|{x}, next(iter(rest)))
        if c is None or len(c)!=len(rest): return False
    return True

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 2311)
byD=collections.defaultdict(lambda:[0,0]); byQ=collections.defaultdict(lambda:[0,0])
nb=0
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 4):
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    nb+=1
    Ps=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        Ps += [(i,j,p) for p in rand_paths(a,terms[i],terms[j],rng,1200) if (len(p)-1)%2==0]
    Ds=[]; seen=set()
    for _ in range(1200):
        c=rand_cycle(a,rng)
        if c and len(c)%2==0 and frozenset(c) not in seen:
            seen.add(frozenset(c)); Ds.append(c)
    Ds.sort(key=len)
    for (i,j,P) in Ps[:100]:
        for D in Ds[:100]:
            S=set(P)&set(D)
            if len(S)<2: continue
            run=arc_span(D,S)
            if run is None: continue
            u,w=run[0],run[-1]
            if P.index(u)>P.index(w): u,w=w,u
            n=len(D); iu=D.index(u)
            for cand in (D[(iu-1)%n], D[(iu+1)%n]):
                if cand in S: continue
                ys=[y for y in a[cand] if y not in D]
                if not ys: continue
                yp=ys[0]; Q2=set(P[P.index(w):])
                if cand in Q2: continue
                drop=set(D)|Q2
                if terms[i] in drop or yp in drop: continue
                C=comp_of(a,drop,terms[i])
                good = bool(C and yp in C and nonbip(a,C) and two_conn(a,C))
                byD[len(D)][good]+=1; byQ[min(len(Q2),14)][good]+=1
                break
print('  bridges: %d'%nb)
print('  |D|    fires / total    rate')
for k in sorted(byD):
    f,t=byD[k][1],byD[k][0]+byD[k][1]
    print('   %2d     %5d / %-6d   %.3f'%(k,f,t,f/max(t,1)))
print('  |Q2|   fires / total    rate')
for k in sorted(byQ):
    f,t=byQ[k][1],byQ[k][0]+byQ[k][1]
    print('   %2d     %5d / %-6d   %.3f'%(k,f,t,f/max(t,1)))
