"""minusD342.py -- (PT2) reduces to: B - D non-bipartite and connected.

To shift a clean meeting's arc endpoint from u to the adjacent u' on D, the path must reach
u' by its unique third edge, from y' say, and the piece from b_i to y' must AVOID D.  Lemma
B (326.2) supplies that piece in EITHER parity provided B - D is connected and
non-bipartite; choosing the parity that keeps P even gives an arc of length |A| +/- 1, of
opposite parity, hence (PT2).
So: is B - D non-bipartite and connected, for the even cycles D that admit clean meetings?
Note B - D non-bipartite <=> some ODD CYCLE IS DISJOINT FROM D, which ties this to the
odd-cycle packing strand of 223 and 254.
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
def sub(a, drop):
    b={}
    for u in a:
        if u in drop: continue
        b[u]={w for w in a[u] if w not in drop}
    return b
def bip(a):
    col={}
    for s in a:
        if s in col: continue
        col[s]=0;q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in col: col[w]=1-col[u];q.append(w)
                elif col[w]==col[u]: return False
    return True
def conn(a):
    if not a: return False
    s=next(iter(a)); seen={s}; q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w not in seen: seen.add(w); q.append(w)
    return len(seen)==len(a)
def arc_len(D,S):
    n=len(D); idx=[i for i,x in enumerate(D) if x in S]
    if len(idx)!=len(S) or not idx: return None
    for s in idx:
        if set(D[(s+k)%n] for k in range(len(S)))==S: return len(S)-1
    return None

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 1777)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 5):
    a,terms=make_real_bridge(30,rng)
    if a is None: st['generation failed']+=1; continue
    st['genuine bridges (n=%d)'%len(a)]+=1
    Ps=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        Ps += [p for p in rand_paths(a,terms[i],terms[j],rng,2000) if (len(p)-1)%2==0]
    seen=set(); good=0; tot=0
    for _ in range(2000):
        c=rand_cycle(a,rng)
        if not c or len(c)%2: continue
        k=frozenset(c)
        if k in seen: continue
        seen.add(k)
        clean=any(arc_len(c,set(P)&set(c)) is not None and len(set(P)&set(c))>=2
                  for P in Ps[:150])
        if not clean: continue
        tot+=1
        R=sub(a,set(c)-set(terms))
        R={u:ws for u,ws in R.items()}
        ok = conn(R) and not bip(R)
        if ok: good+=1
    st['  even cycles with a clean meeting: %d'%tot]+=1
    st['  of which B-D connected AND non-bipartite: %d'%good]+=1
    st['CONDITION HOLDS for some D' if good else '*** condition fails for all D ***']+=1
for k in sorted(st): print('  %-52s %s'%(k,st[k]))
