"""pos375.py -- measure (BP2): the checkerboard refutation, exactly per attachment.

POSITION LEMMA: (SC) + even connecting path P (pair {i,j}) + path Q from b_k attaching
at INTERIOR z of P (internally disjoint from P) with d(z)+|Q| EVEN forces
    2 d(z) == c_ij + c_ik - c_jk (mod 4),  i.e.  d(z) mod 2 is a GLOBAL constant.
A vertex z with b_k-paths of BOTH parities in B - (P - z) realises an even-completing
attachment (and a balanced fan).  So TWO both-parity attachment vertices at OPPOSITE
position parities on one P refute (SC):  (BP2).
Measure: for sampled even connecting paths, the both-parity attachment vertices and
their position parities -- availability exact via Lemma B' block-chains.
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
    if s not in verts or t not in verts: return (False,False,None)
    comp={s}; q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w in verts and w not in comp: comp.add(w); q.append(w)
    if t not in comp: return (False,False,None)
    if s==t: return (True,False,0)
    blks=[b for b in blocks_of(a,comp)]
    incid=collections.defaultdict(list)
    for bi,b in enumerate(blks):
        for v in b: incid[v].append(bi)
    cuts={v for v in incid if len(incid[v])>1}
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

def rand_conn(a, terms, i, j, rng):
    """random b_i - b_j path avoiding the third terminal internally."""
    k=3-i-j; s,t=terms[i],terms[j]; bad=terms[k]
    p=[s]; on={s,bad}
    while p[-1]!=t:
        nb=[w for w in a[p[-1]] if w not in on]
        if not nb: return None
        w=rng.choice(nb); p.append(w); on.add(w)
        if len(p)>len(a): return None
    return p

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 375)
NB=int(sys.argv[2]) if len(sys.argv)>2 else 8
PATHS=int(sys.argv[3]) if len(sys.argv)>3 else 40
made=0
print('  per bridge: even paths tested / (BP2) witnesses / both-parity z per path (min-max)')
while made<NB:
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    if not tau_ge2(a): continue
    made+=1
    wit=0; tested=0; bpcounts=[]; parity_profiles=collections.Counter()
    for i,j in [(0,1),(0,2),(1,2)]:
        k=3-i-j; got=set()
        for _ in range(3000):
            if len(got)>=PATHS//3+1: break
            P=rand_conn(a,terms,i,j,rng)
            if P is None or (len(P)-1)%2: continue
            key=tuple(P)
            if key in got or tuple(reversed(P)) in got: continue
            got.add(key); tested+=1
            verts0=set(a)-set(P)
            bp=[]
            for d in range(1,len(P)-1):
                z=P[d]
                verts=verts0|{z}
                rc,both,delta=parity_set(a,verts,z,terms[k])
                if rc and both: bp.append(d%2)
            bpcounts.append(len(bp))
            prof=frozenset(bp)
            parity_profiles[tuple(sorted(prof))]+=1
            if prof=={0,1}: wit+=1
    print('  n=33  paths %3d   BP2 %3d   both-parity z per path %d-%d   profiles %s'
          %(tested,wit,min(bpcounts) if bpcounts else -1,
            max(bpcounts) if bpcounts else -1,dict(parity_profiles)))
