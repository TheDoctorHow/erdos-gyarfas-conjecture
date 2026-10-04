"""chain374.py -- EXACT variant status for the Trim Lemma, via Lemma B' (342.2).

For balanced fan F, variant (i,j): X = (L_i u L_j) - w.  The trim fires iff b_k is
reachable from w in B - X and sigma+1 lies in the parity set Pi of w-to-b_k PATHS.
By Lemma B': Pi = {0,1} iff some block on the w-b_k block-chain of B - X is
non-bipartite ("chain-odd"); else Pi = {delta}, delta = sum of 2-colour differences
across the bipartite chain blocks.  All exact -- no randomized R search.
Classes per variant:
   chain-odd   : fires via a non-bipartite chain block         (the (CHAIN) witness)
   lucky       : all chain blocks bipartite but delta == sigma+1  (fires anyway)
   locked      : all chain blocks bipartite, delta == sigma       (cannot fire)
   unreached / short / no-pivot
Fans randomized (positive witnesses only); per-fan analysis exact.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0])

def bip_col(a, verts):
    """2-colouring of induced subgraph, or None if non-bipartite (per component)."""
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

def blocks_of(a, verts):
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

def parity_set(a, verts, s, t):
    """(reachable, both, delta): exact parity set of s-t paths in induced subgraph."""
    if s not in verts or t not in verts: return (False,False,None)
    comp={s}; q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w in verts and w not in comp: comp.add(w); q.append(w)
    if t not in comp: return (False,False,None)
    if s==t: return (True,False,0)
    blks=[b for b in blocks_of(a,comp)]
    # block-cut tree
    incid=collections.defaultdict(list)
    for bi,b in enumerate(blks):
        for v in b: incid[v].append(bi)
    cuts={v for v in incid if len(incid[v])>1}
    # tree nodes: ('B',i) and ('C',v); BFS from a block containing s to one containing t
    def nodes_of(v): return [('B',i) for i in incid[v]]
    start=[('B',i) for i in incid[s]]; goal=set(('B',i) for i in incid[t])
    prev={}; q=deque(start)
    for nd in start: prev[nd]=None
    endnode=None
    while q:
        nd=q.popleft()
        if nd in goal: endnode=nd; break
        if nd[0]=='B':
            for v in blks[nd[1]]:
                if v in cuts and ('C',v) not in prev:
                    prev[('C',v)]=nd; q.append(('C',v))
        else:
            for bi in incid[nd[1]]:
                if ('B',bi) not in prev:
                    prev[('B',bi)]=nd; q.append(('B',bi))
    chain=[]; nd=endnode
    while nd is not None: chain.append(nd); nd=prev[nd]
    chain.reverse()
    both=False; delta=0; entry=s
    for pos,nd in enumerate(chain):
        if nd[0]!='B': continue
        blk=set(blks[nd[1]])
        exit_v = t if pos==len(chain)-1 else chain[pos+1][1]
        if entry==exit_v: continue
        col=bip_col(a,blk)
        if col is None: both=True
        else: delta^= (col[entry]!=col[exit_v])
        entry=exit_v
    return (True,both,None if both else delta)

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 374)
NB=int(sys.argv[2]) if len(sys.argv)>2 else 8
FANS=int(sys.argv[3]) if len(sys.argv)>3 else 60
made=0
print('  per bridge: fans / chain-odd fans / firing fans (any route) / variant classes')
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
    codd=0; fire=0; tal=collections.Counter()
    for v,legs,mus in fans:
        sigma=mus[0]%2; hit_odd=False; hit=False
        for i,j in itertools.permutations((0,1,2),2):
            k=3-i-j
            if mus[i]<2: tal['short']+=1; continue
            w=legs[i][1]
            if all(y in legs[i] for y in a[w]): tal['no-pivot']+=1; continue
            X=(set(legs[i])|set(legs[j]))-{w}
            verts=set(a)-X
            rc,both,delta=parity_set(a,verts,w,terms[k])
            if not rc: tal['unreached']+=1
            elif both: tal['chain-odd']+=1; hit_odd=True; hit=True
            elif delta==(sigma+1)%2: tal['lucky']+=1; hit=True
            else: tal['locked']+=1
        if hit_odd: codd+=1
        if hit: fire+=1
    print('  n=%2d  fans %3d   chain-odd fans %3d   firing fans %3d   %s'
          %(len(a),len(fans),codd,fire,dict(tal)))
