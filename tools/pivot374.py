"""pivot374.py -- measure (PIVOT-ODD), the Detour Corollary's premise, exactly.

For balanced fan F, variant (i,j): X = (L_i union L_j) - w, x = pivot's third neighbour.
  (a) b_k reachable from w in B - X                       [BFS]
  (b) x lies on an odd cycle of H := B - (L_i u L_j u w)  [exact: some block of H
      containing x is non-bipartite, by Lemma A inside the block]
(PIVOT-ODD): some fan + variant has (a) and (b).  373.3: (PIVOT-ODD) ==> trim fires
==> not (SC).  Fan sampling randomized (positive witnesses only); (a),(b) exact per fan.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0])

def bipartite_sub(a, verts):
    col={}
    for s in verts:
        if s in col: continue
        col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in verts: continue
                if w not in col: col[w]=1-col[u]; q.append(w)
                elif col[w]==col[u]: return False
    return True

def tau_ge2(a):
    if not bipartite_sub(a,set(a)): pass
    else: return False
    for e in {frozenset((u,w)) for u in a for w in a[u]}:
        u,v=tuple(e)
        b={x:set(a[x]) for x in a}; b[u].discard(v); b[v].discard(u)
        if bipartite_sub(b,set(b)): return False
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

def reach(a, drop, s, t):
    if s in drop or t in drop: return False
    seen={s}; q=deque([s])
    while q:
        u=q.popleft()
        if u==t: return True
        for w in a[u]:
            if w not in drop and w not in seen: seen.add(w); q.append(w)
    return False

def blocks_of(a, verts):
    """biconnected components (as vertex-frozensets) of the induced subgraph on verts."""
    adj={u:[w for w in a[u] if w in verts] for u in verts}
    seen=set(); out=[]
    for root in verts:
        if root in seen: continue
        num={}; low={}; stack=[(root,None,iter(adj[root]))]; num[root]=low[root]=0
        cnt=1; estack=[]
        while stack:
            u,pu,it=stack[-1]
            adv=False
            for w in it:
                if w==pu: continue
                if w not in num:
                    estack.append((u,w)); num[w]=low[w]=cnt; cnt+=1
                    stack.append((w,u,iter(adj[w]))); adv=True; break
                elif num[w]<num[u]:
                    estack.append((u,w)); low[u]=min(low[u],num[w])
            if adv: continue
            stack.pop()
            if stack:
                p=stack[-1][0]
                low[p]=min(low[p],low[u])
                if low[u]>=num[p]:
                    comp=set()
                    while estack:
                        e=estack.pop(); comp.add(e[0]); comp.add(e[1])
                        if e==(p,u): break
                    if comp: out.append(frozenset(comp))
        seen|=set(num)
    return out

def odd_cycle_through(a, verts, x):
    """exact: x on an odd cycle of induced subgraph on verts."""
    if x not in verts: return False
    for blk in blocks_of(a, verts):
        if x in blk and len(blk)>=3 and not bipartite_sub(a, set(blk)):
            return True
    return False

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 374)
NB=int(sys.argv[2]) if len(sys.argv)>2 else 8
FANS=int(sys.argv[3]) if len(sys.argv)>3 else 60
made=0
print('  per bridge: fans / fans with (a)&(b) in some variant / variant tallies')
while made<NB:
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    if not tau_ge2(a): continue
    made+=1
    interior=[v for v in a if v not in terms and len(a[v])==3]
    fans=[]; seenf=set()
    for _ in range(4000):
        if len(fans)>=FANS: break
        v=rng.choice(interior)
        legs=sample_fan_legs(a,terms,v,rng)
        if legs is None: continue
        mus=[len(legs[t])-1 for t in (0,1,2)]
        if not (mus[0]%2==mus[1]%2==mus[2]%2): continue
        key=(v,tuple(tuple(legs[t]) for t in (0,1,2)))
        if key in seenf: continue
        seenf.add(key); fans.append((v,legs,mus))
    ok=0; tal=collections.Counter()
    for v,legs,mus in fans:
        hit=False
        for i,j in itertools.permutations((0,1,2),2):
            k=3-i-j
            if mus[i]<2: tal['short']+=1; continue
            w=legs[i][1]
            x=[y for y in a[w] if y!=v and y!=legs[i][2 if mus[i]>=2 else 1]]
            x=[y for y in a[w] if y not in legs[i]]
            if not x: tal['no-pivot']+=1; continue
            x=x[0]
            X=(set(legs[i])|set(legs[j]))-{w}
            ra=reach(a,X,w,terms[k])
            H=set(a)-set(legs[i])-set(legs[j])-{w}
            rb=odd_cycle_through(a,H,x)
            tal[('a' if ra else '-')+('b' if rb else '-')]+=1
            if ra and rb: hit=True
        if hit: ok+=1
    print('  n=%2d  fans %3d   (PIVOT-ODD) fans %3d   %s'
          %(len(a),len(fans),ok,dict(tal)))
