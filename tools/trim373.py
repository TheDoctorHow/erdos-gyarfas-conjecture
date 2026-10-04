"""trim373.py -- measure (TRIM-EX): does some balanced fan admit a trim?

TRIM LEMMA (arithmetic): F balanced fan at v, legs L_1,L_2,L_3, parity sigma, mu_1 >= 2,
w = first vertex of L_1.  If a path R from w to b_3 exists with |R| == sigma+1 (mod 2),
internally disjoint from tail(L_1) and from {v} union L_2, then
    P_13 = tail(L_1)+R  and  P_23 = L_2+(v,w)+R  are EVEN connecting paths, and (SC)
    demands q == mu_3 + 1  AND  q == mu_3 - 1 (mod 4) for q = |R|: contradiction.
Six variants per fan (ordered pairs (i,j), k the third).  Measure: fraction of balanced
fans admitting >= 1 trim; whether every bridge has some trim-admitting fan; failure mode
split (wrong parity only vs unreachable).  Randomized: positive witnesses only.
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

def sample_fan_legs(a, terms, v, rng):
    tset=set(terms); blocked={v}; legs={}
    for t in rng.sample([0,1,2],3):
        tgt=terms[t]; p=[v]; on=set(blocked)
        while p[-1]!=tgt:
            nb=[w for w in a[p[-1]] if w not in on and (w==tgt or w not in tset)]
            if not nb: return None
            w=rng.choice(nb); p.append(w); on.add(w)
            if len(p)>len(a): return None
        legs[t]=p; blocked|=set(p)-{v}
    return legs

def trim_parities(a, w, target, avoid, rng, K):
    """parities of w->target path lengths found, interior avoiding 'avoid'."""
    found=set()
    for _ in range(K):
        p=[w]; on={w}|avoid
        while p[-1]!=target:
            nb=[x for x in a[p[-1]] if x not in on]
            if not nb: break
            x=rng.choice(nb); p.append(x); on.add(x)
            if len(p)>len(a): break
        if p[-1]==target: found.add((len(p)-1)%2)
        if found=={0,1}: break
    return found

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 373)
NB=int(sys.argv[2]) if len(sys.argv)>2 else 8
FANS=int(sys.argv[3]) if len(sys.argv)>3 else 60
K=int(sys.argv[4]) if len(sys.argv)>4 else 250
made=0
print('  per bridge: balanced fans tested / admitting a trim / variant outcomes')
while made<NB:
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    if not tau_ge2(a): continue
    made+=1
    interior=[v for v in a if v not in terms and len(a[v])==3]
    fans=[]; seen=set()
    for _ in range(4000):
        if len(fans)>=FANS: break
        v=rng.choice(interior)
        legs=sample_fan_legs(a,terms,v,rng)
        if legs is None: continue
        mus=[len(legs[t])-1 for t in (0,1,2)]
        if not (mus[0]%2==mus[1]%2==mus[2]%2): continue
        key=(v,tuple(mus),tuple(tuple(legs[t]) for t in (0,1,2)))
        if key in seen: continue
        seen.add(key); fans.append((v,legs,mus))
    ok=0; outc=collections.Counter()
    for v,legs,mus in fans:
        sigma=mus[0]%2; hit=False
        for i,j in itertools.permutations((0,1,2),2):
            k=3-i-j
            if mus[i]<2: outc['short']+=1; continue
            w=legs[i][1]
            avoid=(set(legs[i])|set(legs[j]))-{w}
            ps=trim_parities(a,w,terms[k],avoid,rng,K)
            if (sigma+1)%2 in ps: outc['fires']+=1; hit=True
            elif ps: outc['wrong-parity']+=1
            else: outc['unreached']+=1
        if hit: ok+=1
    print('  n=%2d  fans %3d   trim-admitting %3d   %s'
          %(len(a),len(fans),ok,dict(outc)))
