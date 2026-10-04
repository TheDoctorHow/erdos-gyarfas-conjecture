"""seed378.py -- check the SEED COROLLARY: every bridge with tau >= 2 has an odd cycle
through AT MOST ONE terminal.  (Proved in 378 via the death of W2; this is the hygiene
check.)  Exact: BFS-parity per terminal-deleted subgraph.
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

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 378)
made=0; ok0=0; ok1=0
while made<12:
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    if not tau_ge2(a): continue
    made+=1
    # odd cycle avoiding ALL terminals?
    free = bip_col(a, set(a)-set(terms)) is None
    # odd cycle through exactly one terminal b_m = odd cycle in B - other two
    one=False
    for m in range(3):
        oth=[terms[x] for x in range(3) if x!=m]
        if bip_col(a, set(a)-set(oth)) is None and not (bip_col(a,set(a)-set(terms)) is None and False):
            # non-bipartite after deleting the other two; if also terminal-free odd exists
            # this still witnesses <=1; for EXACTLY one we would subtract, but <=1 is the claim
            one=True
    if free: ok0+=1
    if free or one: ok1+=1
print('bridges tau>=2 tested: %d   terminal-free odd cycle: %d   odd cycle meeting <=1 terminal: %d'
      %(made,ok0,ok1))
