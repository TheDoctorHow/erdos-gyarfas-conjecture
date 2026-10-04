"""lock372.py -- measure support for (LOCK-BREAK).

(SC) pins the full invariant  I(F) = ((mu1-mu2) mod 4, (mu1-mu3) mod 4, mu1 mod 2)
of every balanced fan F to constants of B:
    mu1-mu2 == c13-c23,  mu1-mu3 == c12-c23  (mod 4),  2*mu1 == c12+c13-c23 (mod 4).
So TWO balanced fans with different I refute (SC).  (BAL) targeted only the third
coordinate (needs opposite parities; the even side binds).  Here: do real bridges with
tau_e >= 2 carry balanced fans with >= 2 distinct I-values, and in particular among
SAME-parity (odd) fans only, and at a SINGLE centre?
Random fan sampling: positive existence only; misses prove nothing.  Self-contained.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0])

def bipartite(a, skip_edge=None):
    col={}
    for s in a:
        if s in col: continue
        col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if skip_edge and frozenset((u,w))==skip_edge: continue
                if w not in col: col[w]=1-col[u]; q.append(w)
                elif col[w]==col[u]: return False
    return True

def tau_ge2(a):
    if bipartite(a): return False
    E={frozenset((u,w)) for u in a for w in a[u]}
    for e in E:
        if bipartite(a, skip_edge=e): return False
    return True

def sample_fan(a, terms, v, rng):
    """three internally disjoint paths v -> b0,b1,b2, terminals internal to no leg.
    Returns (mu0,mu1,mu2) or None."""
    tset=set(terms); blocked={v}; mu={}
    for t in rng.sample([0,1,2],3):
        tgt=terms[t]; p=[v]; on=set(blocked)
        while p[-1]!=tgt:
            nb=[w for w in a[p[-1]] if w not in on and (w==tgt or w not in tset)]
            if not nb: return None
            w=rng.choice(nb); p.append(w); on.add(w)
            if len(p)>len(a): return None
        mu[t]=len(p)-1; blocked|=set(p)-{v}
    return (mu[0],mu[1],mu[2])

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 372)
NB=int(sys.argv[2]) if len(sys.argv)>2 else 8
TRIES=int(sys.argv[3]) if len(sys.argv)>3 else 300
made=0
print('  per bridge: balanced fans found / distinct I overall / distinct I among sigma=1 /')
print('  among sigma=0 / max distinct I at one centre / sigmas present')
while made<NB:
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    if not tau_ge2(a): continue
    made+=1
    interior=[v for v in a if v not in terms and len(a[v])==3]
    seen=set()                      # (v, mu-triple) dedup
    Iall=set(); I1=set(); I0=set()
    bycen=collections.defaultdict(set)
    nfan=0
    for v in interior:
        for _ in range(TRIES):
            f=sample_fan(a,terms,v,rng)
            if f is None: continue
            m0,m1,m2=f
            if not (m0%2==m1%2==m2%2): continue
            if (v,f) in seen: continue
            seen.add((v,f)); nfan+=1
            I=((m0-m1)%4,(m0-m2)%4,m0%2)
            Iall.add(I); bycen[v].add(I)
            (I1 if m0%2 else I0).add(I)
    maxc=max((len(s) for s in bycen.values()), default=0)
    sig=sorted({i[2] for i in Iall})
    print('  n=%2d  fans %5d   |I| %2d   |I,s=1| %2d   |I,s=0| %2d   maxcentre %2d   sigmas %s'
          %(len(a),nfan,len(Iall),len(I1),len(I0),maxc,sig))
