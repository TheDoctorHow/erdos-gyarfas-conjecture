"""diag318.py -- clean diagnostic at FULL search strength: is |N| the governing parameter?
For each bridge with tau>=2 use the MINIMUM-|N| colouring with monochromatic terminals,
count the vertices carrying a uniform-parity Y, and record whether two are adjacent.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/yparity315.py').read().split('rng=random.Random')[0].split('"""',2)[2])

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 41)
st=collections.Counter(); rows=[]
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 40):
    a,terms=make_bridge(rng.choice([10,12,14]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    t=tau_edge(E)
    if t<2: continue
    best=None
    for trial in range(250):
        col={v:rng.randint(0,1) for v in a}
        for x in terms: col[x]=0
        imp=True
        while imp:
            imp=False
            for v in a:
                if v in terms: continue
                nc=sum(1 for w in a[v] if col[w]==col[v])
                if nc>len(a[v])-nc: col[v]^=1; imp=True
        N={e for e in E if col[tuple(e)[0]]==col[tuple(e)[1]]}
        if best is None or len(N)<len(best[1]): best=(dict(col),set(N))
    col,N=best
    P={v:y_parities(a,v,terms,N,rng,trials=160) for v in a if v not in terms}
    uni={v for v in P if (0,0,0) in P[v] or (1,1,1) in P[v]}
    hit=any((not (e in N)) and (tuple(e)[0] in uni) and (tuple(e)[1] in uni)
            for e in E if not (set(tuple(e)) & set(terms)))
    st['bridges with tau>=2']+=1
    st['  |N|=%d : %s'%(len(N),'FIRES' if hit else 'miss')]+=1
    rows.append((len(N),len(uni),len(a)-3,hit))
for k in sorted(st): print('  %-34s %s'%(k,st[k]))
byN=collections.defaultdict(list)
for n_,u_,tot,h in rows: byN[n_].append((u_/max(tot,1),h))
print('  |N|   mean fraction of vertices with uniform-parity Y   fire rate')
for n_ in sorted(byN):
    v=byN[n_]
    print('   %2d        %.2f                                    %d/%d'
          %(n_, sum(x for x,_ in v)/len(v), sum(1 for _,h in v if h), len(v)))
