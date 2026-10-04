"""gencheck338.py -- do the generated bridges satisfy the programme's own hypotheses?

213.1 requires the three terminals to be PAIRWISE AT DISTANCE AT LEAST 5.  The (CYC)
witnesses of 336 had arcs of length 2 and 4, which would put two terminals at distance 2.
Check directly: terminal degrees, pairwise distances, girth, subcubicity.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split(chr(34)*3,2)[2])

def dist(a,s):
    d={s:0};q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w not in d: d[w]=d[u]+1;q.append(w)
    return d
def girth(a):
    best=10**9
    for s in a:
        d={s:0};par={s:None};q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in d: d[w]=d[u]+1;par[w]=u;q.append(w)
                elif w!=par[u]: best=min(best,d[u]+d[w]+1)
    return best

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 401)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 40):
    a,terms=make_bridge(rng.choice([10,12,14]),rng)
    if a is None: continue
    st['bridges generated']+=1
    D=[dist(a,t) for t in terms]
    pd=[D[i][terms[j]] for i,j in [(0,1),(0,2),(1,2)]]
    st['  terminal pairwise distances %s'%str(tuple(sorted(pd)))]+=1
    st['  min terminal distance = %d'%min(pd)]+=1
    st['SATISFIES 213.1 (all pairwise >= 5)' if min(pd)>=5
       else '*** VIOLATES 213.1 ***']+=1
    st['  girth = %d'%girth(a)]+=1
    st['  terminal degrees %s'%str(tuple(sorted(len(a[t]) for t in terms)))]+=1
for k in sorted(st): print('  %-46s %s'%(k,st[k]))
