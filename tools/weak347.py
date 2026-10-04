"""weak347.py -- the weakened Hexagon Theorem hypothesis, measured DIRECTLY.

The proof of 344.2 uses only:
  (i)  D even, through a terminal b_j  (so the tail is empty);
  (ii) TWO ADJACENT vertices x, x' of D whose third edges leave D  (weaker than chordless);
  (iii) b_i-to-y and b_i-to-y' paths of BOTH PARITIES inside B - D.
By Lemma B' (342.2), (iii) is EXACTLY "some block on the chain is non-bipartite", so rather
than test a proxy - the error of 331 and 341 - test (iii) DIRECTLY by exhibiting paths of
both parities.  A positive find is conclusive.
Also splits the failures: disconnection versus single-parity.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def both_parities(a, drop, s, t, rng, K=400):
    par=set()
    for _ in range(K):
        p=[s]; on={s}|set(drop)
        if s in drop or t in drop: return None
        while p[-1]!=t:
            nb=[w for w in a[p[-1]] if w not in on]
            if not nb: break
            w=rng.choice(nb); p.append(w); on.add(w)
        if p[-1]==t:
            par.add((len(p)-1)%2)
            if len(par)==2: return True
    return False if par else None

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 3121)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 5):
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    st['genuine bridges']+=1
    seen=set(); cand=[]
    for _ in range(5000):
        c=rand_cycle(a,rng)
        if not c or len(c)%2: continue
        k=frozenset(c)
        if k in seen: continue
        seen.add(k)
        if any(t in c for t in terms): cand.append(c)
    ok=0; noreach=0; onepar=0; noadj=0
    for D in cand:
        Ds=set(D); n=len(D)
        bj=[t for t in terms if t in Ds][0]
        bi=[t for t in terms if t not in Ds]
        if not bi: continue
        bi=bi[0]
        pair=None
        for k in range(n):
            x,xp=D[k],D[(k+1)%n]
            ys=[y for y in a[x] if y not in Ds]; yps=[y for y in a[xp] if y not in Ds]
            if ys and yps: pair=(ys[0],yps[0]); break
        if pair is None: noadj+=1; continue
        r1=both_parities(a,Ds,bi,pair[0],rng)
        r2=both_parities(a,Ds,bi,pair[1],rng)
        if r1 is None or r2 is None: noreach+=1
        elif r1 and r2: ok+=1
        else: onepar+=1
    st['  terminal even cycles: %d'%len(cand)]+=1
    st['  HYPOTHESIS (iii) VERIFIED directly: %d'%ok]+=1
    st['  failed: no adjacent pair with outside edges: %d'%noadj]+=1
    st['  failed: y unreachable in B-D: %d'%noreach]+=1
    st['  failed: only one parity available: %d'%onepar]+=1
    st['THEOREM APPLIES' if ok else '*** no qualifying cycle ***']+=1
for k in sorted(st): print('  %-54s %s'%(k,st[k]))
