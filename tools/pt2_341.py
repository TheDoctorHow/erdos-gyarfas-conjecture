"""pt2_341.py -- (PT2): one even cycle, two clean meetings, arcs of different parity.

The parity test gives |A| == m (mod 2) for every clean meeting of an even cycle D = 2m by
an even connecting path.  So if a SINGLE even cycle admits two clean meetings whose arcs
differ in parity, (SC) fails outright - and |D| mod 4 need never be computed.
     (PT2)  some even cycle admits clean meetings with arcs of BOTH parities.
Strictly weaker still than (PT): no arithmetic on |D|, only a parity clash between two arcs.
Tested on genuine bridges (213.1 verified, n = 33).
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
def arc_len(D,S):
    n=len(D); idx=[i for i,x in enumerate(D) if x in S]
    if len(idx)!=len(S) or not idx: return None
    for s in idx:
        if set(D[(s+k)%n] for k in range(len(S)))==S: return len(S)-1
    return None

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 1511)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 5):
    a,terms=make_real_bridge(30,rng)
    if a is None: st['generation failed']+=1; continue
    st['genuine bridges (n=%d)'%len(a)]+=1
    Ps=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        Ps += [p for p in rand_paths(a,terms[i],terms[j],rng,3000) if (len(p)-1)%2==0]
    Ds=[]
    seen=set()
    for _ in range(3000):
        c=rand_cycle(a,rng)
        if c and len(c)%2==0:
            k=frozenset(c)
            if k not in seen: seen.add(k); Ds.append(c)
    hit=None
    for D in Ds[:250]:
        par=collections.defaultdict(list)
        for P in Ps[:250]:
            S=set(P)&set(D)
            if len(S)<2: continue
            L=arc_len(D,S)
            if L is None: continue
            par[L%2].append(L)
        if 0 in par and 1 in par:
            hit=(len(D),sorted(set(par[0]))[:3],sorted(set(par[1]))[:3]); break
    st['*** (PT2) WITNESS -> (SC) fails ***' if hit else 'no parity clash found (inconclusive)']+=1
    if hit: st['  |D|=%d  even arcs %s  odd arcs %s'%hit]+=1
    st['  even cycles sampled %d, even connecting paths %d'%(len(Ds),len(Ps))]+=1
for k in sorted(st): print('  %-54s %s'%(k,st[k]))
