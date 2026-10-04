"""full348.py -- (TERM) with the freedom it actually has.

345.2 fixed b_i as the first terminal off D and took the first adjacent pair (x,x') with
outside third edges.  (TERM) allows BOTH choices of b_i and ALL such adjacent pairs.  A
cycle qualifies if ANY combination gives parity-flexible reachability to both outside
neighbours.  Measured on genuine bridges, hypothesis tested directly.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def both_parities(a, drop, s, t, rng, K=300):
    if s in drop or t in drop: return None
    par=set()
    for _ in range(K):
        p=[s]; on={s}|set(drop)
        while p[-1]!=t:
            nb=[w for w in a[p[-1]] if w not in on]
            if not nb: break
            w=rng.choice(nb); p.append(w); on.add(w)
        if p[-1]==t:
            par.add((len(p)-1)%2)
            if len(par)==2: return True
    return False if par else None

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 3313)
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
    first=0; full=0; npairs=[]
    for D in cand:
        Ds=set(D); n=len(D)
        offs=[t for t in terms if t not in Ds]
        pairs=[]
        for k in range(n):
            x,xp=D[k],D[(k+1)%n]
            ys=[y for y in a[x] if y not in Ds]; yps=[y for y in a[xp] if y not in Ds]
            if ys and yps: pairs.append((ys[0],yps[0]))
        npairs.append(len(pairs))
        if not pairs or not offs: continue
        # 345-style: first terminal, first pair
        y,yp=pairs[0]; bi=offs[0]
        if both_parities(a,Ds,bi,y,rng) and both_parities(a,Ds,bi,yp,rng): first+=1
        # full freedom
        got=False
        for bi in offs:
            for (y,yp) in pairs:
                if both_parities(a,Ds,bi,y,rng) and both_parities(a,Ds,bi,yp,rng):
                    got=True; break
            if got: break
        if got: full+=1
    st['  terminal even cycles: %d'%len(cand)]+=1
    st['  qualify, 345-style (first b_i, first pair): %d'%first]+=1
    st['  qualify, FULL freedom: %d'%full]+=1
    st['  mean adjacent pairs available: %.1f'%(sum(npairs)/max(len(npairs),1))]+=1
    st['THEOREM APPLIES' if full else '*** none ***']+=1
for k in sorted(st): print('  %-54s %s'%(k,st[k]))
