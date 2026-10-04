"""pt340.py -- (PT), the parity test, as a target strictly weaker than (CYC).

222's parity test: if an even cycle D is CLEANLY met by an even connecting path P - the
intersection being a single arc A of D - then (SC) forces
        |D| == 2|A|   (mod 4).
Taking D = Z a cycle through two terminals and A = one arc recovers (CYC), so (CYC) is a
SPECIAL CASE and (PT) is the weaker target: one cycle, one path, no two-arcs-even condition
and no second terminal pair.
Tested on GENUINE bridges (213.1 verified, n = 33, above the 30-vertex threshold).
Randomized witness hunt - a find is conclusive.
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
        if v in a[w] and len(p)>=4: return p+[v] if v!=p[0] else p
        if len(p)>len(a): return None

def arc_of(D, S):
    """is the vertex set S a contiguous arc of cycle D (list, D[0]!=D[-1])? return its length"""
    n=len(D); idx=[i for i,x in enumerate(D) if x in S]
    if len(idx)!=len(S) or not idx: return None
    for start in idx:
        run=[D[(start+k)%n] for k in range(len(S))]
        if set(run)==S: return len(S)-1
    return None

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 1201)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 5):
    a,terms=make_real_bridge(30,rng)
    if a is None: st['generation failed']+=1; continue
    st['genuine bridges (n=%d)'%len(a)]+=1
    Ps=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        Ps += [(i,j,p) for p in rand_paths(a,terms[i],terms[j],rng,2500) if (len(p)-1)%2==0]
    Ds=[]
    for _ in range(2500):
        c=rand_cycle(a,rng)
        if c and len(c)%2==0: Ds.append(c)
    found=None; tested=0
    for (i,j,P) in Ps[:400]:
        for D in Ds[:400]:
            S=set(P)&set(D)
            if len(S)<2: continue
            L=arc_of(D,S)
            if L is None: continue
            tested+=1
            if (len(D)-2*L)%4!=0:
                found=(len(D),L,(len(D)-2*L)%4); break
        if found: break
    st['  clean (path,cycle) pairs tested: %d'%tested]+=1
    st['*** (PT) VIOLATED -> (SC) fails ***' if found else 'parity test held on all tested']+=1
    if found: st['  |D|=%d arc=%d  |D|-2|A| = %d mod 4'%found]+=1
for k in sorted(st): print('  %-50s %s'%(k,st[k]))
