"""cert332.py -- beta(v) by enumerating COLOURINGS, not edge subsets.

beta(v) = min number of non-crossing edges over 2-colourings of B' with v and T* the SAME
colour.  That is exactly the least bipartizing set placing them together.  The certificate
of 329.2 forces beta(v) <= 5 for every v, if B'-F_v stays 2-connected.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split(chr(34)*3,2)[2])

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 257)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 12):
    a,terms=make_bridge(rng.choice([8,10]),rng)
    if a is None: continue
    E=[frozenset((u,w)) for u in a for w in a[u]]
    E=list({e for e in E})
    for x in terms: E.append(frozenset(('T*',x)))
    V=sorted(set().union(*[set(e) for e in E]), key=str)
    idx={v:i for i,v in enumerate(V)}
    n=len(V); ti=idx['T*']
    st['bridges']+=1
    beta={v:99 for v in V}
    for mask in range(1<<(n-1)):
        full=(mask<<1) | 0            # fix V[0] colour 0
        col=[(full>>i)&1 for i in range(n)]
        nc=sum(1 for e in E for u,w in [tuple(e)] if col[idx[u]]==col[idx[w]])
        for v in V:
            if col[idx[v]]==col[ti] and nc<beta[v]: beta[v]=nc
    mx=max(beta[v] for v in V if v!='T*')
    for v in V:
        if v!='T*': st['beta(v) = %d'%beta[v]]+=1
    st['max beta over v = %d'%mx]+=1
    st['*** CERTIFICATE REFUTED: beta > 5 ***' if mx>5
       else 'certificate survives here (max beta <= 5)']+=1
for k in sorted(st): print('  %-44s %s'%(k,st[k]))
