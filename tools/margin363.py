"""margin363.py -- the margin in the uncovered case.

Residue: in a bridge with tau_v <= 3, NOT ALL short even terminal cycles are odd-cycle
transversals.  Measure the margin: among terminal HEXAGONS specifically (the shortest, and
the ones 360.1 cannot certify when tau_v = 3), what fraction are transversals, and is any
bridge close to having them all be?
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def bip_after(a,S):
    col={}
    for s in a:
        if s in S or s in col: continue
        col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w in S: continue
                if w not in col: col[w]=1-col[u]; q.append(w)
                elif col[w]==col[u]: return False
    return True
def tau_v_upto(a,cap=3):
    V=sorted(a,key=str)
    for m in range(cap+1):
        for S in itertools.combinations(V,m):
            if bip_after(a,set(S)): return m
    return cap+1
def terminal_hexagons(a,terms):
    out=[]; seen=set()
    for b in terms:
        nb=sorted(a[b],key=str)
        if len(nb)!=2: continue
        p,q=nb
        stack=[(p,[p])]
        while stack:
            u,path=stack.pop()
            if len(path)==5:
                if u==q:
                    c=[b]+path
                    k=frozenset(c)
                    if k not in seen and len(set(c))==6: seen.add(k); out.append(c)
                continue
            for w in a[u]:
                if w==b or w in path: continue
                stack.append((w,path+[w]))
    return out

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 7503)
st=collections.Counter(); rows=[]
low=0
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 90):
    a,terms=make_real_bridge(rng.choice([30,40]),rng)
    if a is None: continue
    tv=tau_v_upto(a,3)
    if tv>3: continue
    low+=1
    H=terminal_hexagons(a,terms)
    trans=[c for c in H if bip_after(a,set(c))]
    rows.append((tv,len(H),len(trans)))
    st['tau_v=%d  hexagons %d, transversal %d  -> gap(B) at 6: %s'
       %(tv,len(H),len(trans),'YES' if len(trans)<len(H) else ('NO HEXAGON' if not H else 'ALL ARE'))]+=1
    if low>=10: break
print('  bridges with tau_v <= 3 examined: %d'%low)
for k in sorted(st): print('   %s  %s'%(k.ljust(58),st[k]))
tot=sum(h for _,h,_ in rows); tr=sum(t for _,_,t in rows)
print('  total terminal hexagons %d, transversal %d (%.0f%%)'%(tot,tr,100*tr/max(tot,1)))
