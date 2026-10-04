"""exact343.py -- the CORRECT hypothesis of the arc-shift construction, measured.

341.2 measured B - D.  The construction of 341.1 needs a b_i-to-y' path avoiding
D UNION Q_2, of a prescribed parity, so what must hold is:
    y' not on D;  u' not on Q_2;
    b_i and y' in the SAME COMPONENT of B - (D union Q_2);
    that component NON-BIPARTITE  (and 2-connected, for Lemma B as proved).
Measured over actual clean meetings on genuine bridges, not over cycles.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/realbridge339.py').read().split('rng=random.Random')[0].split(chr(34)*3,2)[2])

def rand_cycle(a,rng):
    u=rng.choice(list(a)); v=rng.choice(list(a[u]))
    p=[u]; on={u}
    while True:
        nb=[w for w in a[p[-1]] if w not in on and not (len(p)==1 and w==v)]
        if not nb: return None
        w=rng.choice(nb); p.append(w); on.add(w)
        if v in a[w] and len(p)>=4: return p
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
        c=comp_of(a,(set(a)-vs)|{x}, next(iter(vs-{x})))
        if c is None or len(c)!=len(vs)-1: return False
    return True

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 2111)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 4):
    a,terms=make_real_bridge(30,rng)
    if a is None: st['generation failed']+=1; continue
    st['genuine bridges (n=%d)'%len(a)]+=1
    Ps=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        Ps += [(i,j,p) for p in rand_paths(a,terms[i],terms[j],rng,1500) if (len(p)-1)%2==0]
    Ds=[]; seen=set()
    for _ in range(1500):
        c=rand_cycle(a,rng)
        if c and len(c)%2==0 and frozenset(c) not in seen:
            seen.add(frozenset(c)); Ds.append(c)
    tot=0; ok=0; ok2=0
    for (i,j,P) in Ps[:120]:
        for D in Ds[:120]:
            S=set(P)&set(D)
            if len(S)<2: continue
            run=arc_span(D,S)
            if run is None: continue
            u,w=run[0],run[-1]
            if P.index(u)>P.index(w): u,w=w,u
            tot+=1
            n=len(D); iu=D.index(u)
            for cand in (D[(iu-1)%n], D[(iu+1)%n]):
                if cand in S: continue
                up=cand
                ys=[y for y in a[up] if y not in D]
                if not ys: continue
                yp=ys[0]
                Q2=set(P[P.index(w):])
                if up in Q2: continue
                drop=set(D)|Q2
                if terms[i] in drop or yp in drop: continue
                C=comp_of(a,drop,terms[i])
                if C and yp in C and nonbip(a,C):
                    ok+=1
                    if two_conn(a,C): ok2+=1
                break
    st['  clean meetings examined: %d'%tot]+=1
    st['  hypothesis holds (connected+non-bipartite): %d'%ok]+=1
    st['  and component 2-connected: %d'%ok2]+=1
for k in sorted(st): print('  %-54s %s'%(k,st[k]))
