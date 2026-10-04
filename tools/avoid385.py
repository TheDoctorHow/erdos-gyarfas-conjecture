"""avoid385.py -- measure (AVOID) in component form (381.1), exactly per path.

For each pair {i,j} and sampled even connecting path P: is the component of b_k in
B - P NON-BIPARTITE?  (<=> some lollipop avoids P <=> (T2'') witnessed by P.)
Also nu_odd(B) in {1,2,>=3} by greedy disjoint odd cycle search (lower bound only).
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0])

def bip_comp(a, verts, s):
    """is s's component of induced subgraph bipartite? (True/False); None if s absent"""
    if s not in verts: return None
    col={s:0}; q=deque([s]); ok=True
    while q:
        u=q.popleft()
        for w in a[u]:
            if w not in verts: continue
            if w not in col: col[w]=1-col[u]; q.append(w)
            elif col[w]==col[u]: ok=False
    return ok

def tau_ge2(a):
    col={};
    def bip(b):
        c={}
        for s in b:
            if s in c: continue
            c[s]=0; q=deque([s])
            while q:
                u=q.popleft()
                for w in b[u]:
                    if w not in c: c[w]=1-c[u]; q.append(w)
                    elif c[w]==c[u]: return False
        return True
    if bip(a): return False
    for e in {frozenset((u,w)) for u in a for w in a[u]}:
        u,v=tuple(e)
        b={x:set(a[x]) for x in a}; b[u].discard(v); b[v].discard(u)
        if bip(b): return False
    return True

def rand_conn(a, terms, i, j, rng):
    k=3-i-j; s,t=terms[i],terms[j]
    p=[s]; on={s,terms[k]}
    while p[-1]!=t:
        nb=[w for w in a[p[-1]] if w not in on]
        if not nb: return None
        w=rng.choice(nb); p.append(w); on.add(w)
        if len(p)>len(a): return None
    return p

def greedy_nu(a, terms, rng, tries=3000):
    best=0
    for _ in range(30):
        used=set(); cnt=0
        for _ in range(tries):
            c=rand_cycle(a,rng)
            if not c: continue
            cy=c[:-1] if c[0]==c[-1] else c
            if len(cy)%2==0 or set(cy)&used: continue
            used|=set(cy); cnt+=1
        best=max(best,cnt)
    return best

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 385)
made=0
print('  per bridge: nu_odd(greedy) / per pair: even paths tested, (AVOID)-witnesses')
while made<8:
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    if not tau_ge2(a): continue
    made+=1
    nu=greedy_nu(a,terms,rng)
    out=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        k=3-i-j; got=set()
        for _ in range(3000):
            if len(got)>=25: break
            P=rand_conn(a,terms,i,j,rng)
            if P is None or (len(P)-1)%2: continue
            key=tuple(P)
            if key in got or tuple(reversed(P)) in got: continue
            got.add(key)
        wit=sum(1 for key in got if bip_comp(a,set(a)-set(key),terms[k])==False)
        out.append('%d/%d'%(wit,len(got)))
    print('  n=33  nu>=%d   pairs: %s'%(nu,'  '.join(out)))
