"""scale350.py -- do the two gaps of 346.4 persist as the bridge grows?

(i)  How long is the SHORTEST even cycle through a terminal?  346.2 says short cycles always
     qualify, so (TERM-SHORT) needs this to stay small.
(ii) For that shortest cycle, does the Terminal Cycle Theorem's hypothesis hold?
Measured at n = 30, 40, 50 (bridges on n+3 vertices), all above the 30-vertex threshold.
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

def shortest_even_terminal_cycle(a, terms):
    """BFS over (vertex,parity) in B - b for each terminal b with neighbours p,q"""
    best=None
    for b in terms:
        nb=sorted(a[b], key=str)
        if len(nb)!=2: continue
        p,q=nb
        start=(p,0); d={start:0}; Q=deque([start])
        while Q:
            (u,par)=Q.popleft()
            for w in a[u]:
                if w==b: continue
                s=(w,1-par)
                if s not in d: d[s]=d[(u,par)]+1; Q.append(s)
        if (q,0) in d:      # even-length p-q path
            L=d[(q,0)]+2
            if best is None or L<best[0]: best=(L,b,p,q)
    return best

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 3701)
st=collections.Counter()
for NV in (30,40,50):
    for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 3):
        a,terms=make_real_bridge(NV,rng)
        if a is None: st['gen failed n=%d'%NV]+=1; continue
        r=shortest_even_terminal_cycle(a,terms)
        if r is None: st['n=%-3d no even terminal cycle'%len(a)]+=1; continue
        L,b,p,q=r
        st['n=%-3d shortest even terminal cycle = %d'%(len(a),L)]+=1
        # find an actual cycle of that length and test the hypothesis
        got=False
        seen=set()
        for _ in range(8000):
            c=rand_cycle(a,rng)
            if not c or len(c)!=L or b not in c: continue
            k=frozenset(c)
            if k in seen: continue
            seen.add(k)
            Ds=set(c); n=len(c)
            offs=[t for t in terms if t not in Ds]
            for k2 in range(n):
                x,xp=c[k2],c[(k2+1)%n]
                ys=[y for y in a[x] if y not in Ds]; yps=[y for y in a[xp] if y not in Ds]
                if not(ys and yps): continue
                for bi in offs:
                    if both_parities(a,Ds,bi,ys[0],rng) and both_parities(a,Ds,bi,yps[0],rng):
                        got=True; break
                if got: break
            if got or len(seen)>40: break
        st['n=%-3d  theorem applies at shortest cycle: %s'%(len(a),got)]+=1
for k in sorted(st): print('  %s  %s'%(k.ljust(50),st[k]))
