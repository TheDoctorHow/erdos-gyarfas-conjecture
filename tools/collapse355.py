"""collapse355.py -- does the qualification rate depend on |D| or on |D|/n?

349.3 argued the controlling parameter is the FRACTION of the bridge deleted, since what the
Terminal Cycle Theorem needs is B - D to stay reachable and parity-rich.  If so, rates
measured at different n should COLLAPSE when plotted against |D|/n, and separate when plotted
against |D|.  Never tested.  Measured at n = 33 and n = 53.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def both_parities(a,drop,s,t,rng,K=250):
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

def qualifies(a,terms,c,rng):
    Ds=set(c); n=len(c); offs=[t for t in terms if t not in Ds]
    if not offs: return False
    for k in range(n):
        x,xp=c[k],c[(k+1)%n]
        ys=[y for y in a[x] if y not in Ds]; yps=[y for y in a[xp] if y not in Ds]
        if not(ys and yps): continue
        for bi in offs:
            if both_parities(a,Ds,bi,ys[0],rng) and both_parities(a,Ds,bi,yps[0],rng):
                return True
    return False

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 5501)
byL=collections.defaultdict(lambda: collections.defaultdict(lambda:[0,0]))
byR=collections.defaultdict(lambda: collections.defaultdict(lambda:[0,0]))
for NV in (30,50):
    for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 2):
        a,terms=make_real_bridge(NV,rng)
        if a is None: continue
        N=len(a); seen=set()
        for _ in range(4000):
            c=rand_cycle(a,rng)
            if not c or len(c)%2 or not any(t in c for t in terms): continue
            k=frozenset(c)
            if k in seen: continue
            seen.add(k)
            q=qualifies(a,terms,c,rng)
            byL[N][len(c)][1 if q else 0]+=1
            r=round(len(c)/N,1)
            byR[N][r][1 if q else 0]+=1
print('  RATE vs |D|  (absolute length)')
for N in sorted(byL):
    row=' '.join('%d:%.2f'%(L,byL[N][L][1]/max(sum(byL[N][L]),1))
                 for L in sorted(byL[N]) if sum(byL[N][L])>=8)
    print('   n=%-3d %s'%(N,row))
print('  RATE vs |D|/n  (fraction deleted)')
for N in sorted(byR):
    row=' '.join('%.1f:%.2f'%(r,byR[N][r][1]/max(sum(byR[N][r]),1))
                 for r in sorted(byR[N]) if sum(byR[N][r])>=8)
    print('   n=%-3d %s'%(N,row))
