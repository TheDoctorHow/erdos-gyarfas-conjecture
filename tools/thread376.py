"""thread376.py -- measure (THREAD), the theta-threading route to (BP2).

Configuration: odd cycle C avoiding terminals; R_1, R_2 disjoint b_k -> C paths
(except at b_k), internally avoiding C and terminals; T = C u R_1 u R_2; even
connecting path P for the remaining pair {i,j} with P cap T empty; interior vertices
z of P with a path to T - b_k avoiding P and internally avoiding T.
PROVED (376): every s in T - b_k has s-to-b_k paths of BOTH parities inside T, so each
attached z is a both-parity vertex; two attached z at opposite position parities give
(BP2), hence not-(SC).  Greedy-random R construction and P sampling: positive
witnesses only.  Self-contained.
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
    """random-order BFS path from s to any target, interior avoiding banned+targets."""
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

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 376)
NB=int(sys.argv[2]) if len(sys.argv)>2 else 8
made=0
print('  per bridge: (C,k) configs with theta / configs with an even T-free P /')
print('  P tested / (THREAD)-witness P / profile counts')
while made<NB:
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    if not tau_ge2(a): continue
    made+=1
    tset=set(terms)
    # sample distinct odd cycles avoiding terminals
    Cs=[]; seen=set()
    for _ in range(4000):
        c=rand_cycle(a,rng)
        if not c: continue
        cy=c[:-1] if c[0]==c[-1] else c
        if len(cy)%2==0: continue
        if set(cy)&tset: continue
        key=frozenset(cy)
        if key in seen: continue
        seen.add(key); Cs.append(cy)
        if len(Cs)>=6: break
    Cs.sort(key=len)
    cfg_theta=0; cfg_P=0; nP=0; nwit=0; prof=collections.Counter()
    for cy in Cs[:4]:
        Cset=set(cy)
        for k in (0,1,2):
            bk=terms[k]; oth=[terms[m] for m in (0,1,2) if m!=k]
            # two disjoint b_k -> C paths, greedy with retries
            T=None
            for _ in range(40):
                y1,y2=list(a[bk]);
                if rng.random()<0.5: y1,y2=y2,y1
                R1=bfs_path(a,y1,Cset,set(oth)|{bk},rng)
                if R1 is None: continue
                R2=bfs_path(a,y2,Cset,set(oth)|{bk}|set(R1),rng)
                if R2 is None: continue
                T=Cset|set(R1)|set(R2)|{bk}; break
            if T is None: continue
            cfg_theta+=1
            i,j=[m for m in (0,1,2) if m!=k]
            got=set(); found_even=False
            for _ in range(1500):
                if len(got)>=8: break
                P=rand_walk_path(a,terms[i],terms[j],(T|{terms[k]})-{terms[i],terms[j]},rng,len(a))
                if P is None or (len(P)-1)%2: continue
                key=tuple(P)
                if key in got or tuple(reversed(P)) in got: continue
                got.add(key)
            if got: cfg_P+=1
            for key in got:
                P=list(key); Pset=set(P); nP+=1
                allowed=set(a)-Pset-T
                Tgood=T-{bk}
                ppar=set()
                for d in range(1,len(P)-1):
                    z=P[d]
                    # attach: z adjacent to Tgood, or reaches Tgood via allowed
                    hit=any(w in Tgood for w in a[z])
                    if not hit:
                        # BFS in allowed from z's allowed neighbours
                        q=deque([w for w in a[z] if w in allowed]); seenb=set(q)
                        while q and not hit:
                            u=q.popleft()
                            if any(w in Tgood for w in a[u]): hit=True; break
                            for w in a[u]:
                                if w in allowed and w not in seenb:
                                    seenb.add(w); q.append(w)
                    if hit: ppar.add(d%2)
                    if ppar=={0,1}: break
                prof[tuple(sorted(ppar))]+=1
                if ppar=={0,1}: nwit+=1
    print('  n=33  theta-configs %2d   with-P %2d   P %3d   witness %3d   %s'
          %(cfg_theta,cfg_P,nP,nwit,dict(prof)))
