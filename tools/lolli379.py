"""lolli379.py -- (T2)/(T3') under the LOLLIPOP: T = C u R, live set C - {c}.

379: (T1) is PROVED -- Seed Corollary gives an odd cycle C meeting <= 1 terminal; if
one, it is b_k and R is empty (live C - {b_k}); if none, G - e_i - e_j is connected
(3-edge-connectivity), so a b_k -> C path R avoiding b_i, b_j exists (live C - {c}).
Leaves: (T2) an even {i,j}-connecting path P with P cap T empty; (T3') interior
vertices of BOTH position parities of P reaching LIVE = C - {c} off P u T.
Measure exactly per path.  Positive witnesses only for sampled objects.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0])

def bip_col(a, verts):
    col={}
    for s in verts:
        if s in col: continue
        col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in verts: continue
                if w not in col: col[w]=1-col[u]; q.append(w)
                elif col[w]==col[u]: return None
    return col

def tau_ge2(a):
    if bip_col(a,set(a)) is not None: return False
    for e in {frozenset((u,w)) for u in a for w in a[u]}:
        u,v=tuple(e)
        b={x:set(a[x]) for x in a}; b[u].discard(v); b[v].discard(u)
        if bip_col(b,set(b)) is not None: return False
    return True

def bfs_path(a, s, targets, banned, rng):
    if s in banned: return None
    par={s:None}; q=deque([s])
    while q:
        u=q.popleft()
        nbs=list(a[u]); rng.shuffle(nbs)
        for w in nbs:
            if w in par: continue
            if w in targets:
                p=[w,u]
                while par[p[-1]] is not None: p.append(par[p[-1]])
                return list(reversed(p))
            if w in banned: continue
            par[w]=u; q.append(w)
    return None

def rand_walk_path(a, s, t, banned, rng, maxlen):
    p=[s]; on={s}|set(banned)
    while p[-1]!=t:
        nb=[w for w in a[p[-1]] if w==t or w not in on]
        if not nb: return None
        w=rng.choice(nb); p.append(w); on.add(w)
        if len(p)>maxlen: return None
    return p

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 379)
NB=int(sys.argv[2]) if len(sys.argv)>2 else 8
made=0
print('  per bridge: lollipop configs / with even T-free P / P tested / (T3)-witness P / profiles')
while made<NB:
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    if not tau_ge2(a): continue
    made+=1
    tset=set(terms)
    Cs=[]; seen=set()
    for _ in range(4000):
        c=rand_cycle(a,rng)
        if not c: continue
        cy=c[:-1] if c[0]==c[-1] else c
        if len(cy)%2==0 or set(cy)&tset: continue
        key=frozenset(cy)
        if key in seen: continue
        seen.add(key); Cs.append(cy)
        if len(Cs)>=4: break
    Cs.sort(key=len)
    nconf=0; nP_conf=0; nP=0; nwit=0; prof=collections.Counter()
    for cy in Cs[:4]:
        Cset=set(cy)
        for k in (0,1,2):
            bk=terms[k]; oth=[terms[m] for m in (0,1,2) if m!=k]
            R=bfs_path(a,bk,Cset,set(oth),rng)
            if R is None: continue      # cannot happen by 3-edge-connectivity; guard
            c_att=R[-1]
            T=Cset|set(R)
            live=Cset-{c_att}
            nconf+=1
            i,j=[m for m in (0,1,2) if m!=k]
            got=set()
            for _ in range(1500):
                if len(got)>=8: break
                P=rand_walk_path(a,terms[i],terms[j],(T|{bk})-{terms[i],terms[j]},rng,len(a))
                if P is None or (len(P)-1)%2: continue
                key=tuple(P)
                if key in got or tuple(reversed(P)) in got: continue
                got.add(key)
            if got: nP_conf+=1
            for key in got:
                P=list(key); Pset=set(P); nP+=1
                allowed=set(a)-Pset-T
                ppar=set()
                for d in range(1,len(P)-1):
                    z=P[d]
                    hit=any(w in live for w in a[z])
                    if not hit:
                        q=deque([w for w in a[z] if w in allowed]); seenb=set(q)
                        while q and not hit:
                            u=q.popleft()
                            if any(w in live for w in a[u]): hit=True; break
                            for w in a[u]:
                                if w in allowed and w not in seenb:
                                    seenb.add(w); q.append(w)
                    if hit: ppar.add(d%2)
                    if ppar=={0,1}: break
                prof[tuple(sorted(ppar))]+=1
                if ppar=={0,1}: nwit+=1
    print('  n=33  configs %2d   with-P %2d   P %3d   witness %3d   %s'
          %(nconf,nP_conf,nP,nwit,dict(prof)))
