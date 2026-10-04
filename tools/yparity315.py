"""yparity315.py -- (Y-EXIST) weakened to a PARITY condition, and tested.

314.1 required the Y's legs to avoid N.  Redoing the computation without that:
   mu_m == [c(w) =/= c(b_m)] + |L_m cap N|   (mod 2),
so pair (i,j) gives an EVEN connecting path iff |L_i cap N| + c(b_i) == |L_j cap N| + c(b_j).
With the terminals monochromatic that says the three legs have EQUAL N-parity.
And e_m is ODD iff |L_m cap N| == |L'_m cap N| (mod 2), given ww' not in N.
So the mechanism fires as soon as two ADJACENT vertices each carry a Y whose three legs
have the SAME N-parity, and the two Y's agree on that parity.  Y's live in B, not in W.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split('"""', 2)[2])

def bipn(a):
    col={}
    for s in a:
        if s in col: continue
        col[s]=0;q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in col: col[w]=1-col[u];q.append(w)
                elif col[w]==col[u]: return None
    return col
def tau_edge(E, cap=3):
    El=sorted(E,key=lambda s:sorted(map(str,s)))
    for k in range(cap+1):
        for S in itertools.combinations(El,k):
            if bipn(adj_of(set(El)-set(S))) is not None: return k
    return cap+1

def rand_path(a, s, t, banned, rng):
    """random simple s-t path avoiding banned (a set of vertices), or None"""
    d={s:None}; q=deque([s])
    order=lambda xs:(lambda l:(rng.shuffle(l),l)[1])(list(xs))
    while q:
        u=q.popleft()
        if u==t: break
        for x in order(a[u]):
            if x in banned or x in d: continue
            d[x]=u; q.append(x)
    if t not in d: return None
    p=[]; z=t
    while z is not None: p.append(z); z=d[z]
    return p[::-1]

def y_parities(a, w, terms, N, rng, trials=160):
    """achievable vectors (|L_1 cap N|,|L_2 cap N|,|L_3 cap N|) mod 2 for Y's at w"""
    out=set()
    for _ in range(trials):
        for perm in [list(p) for p in itertools.permutations(range(3))]:
            used={w}; legs={}; ok=True
            for i in perm:
                t=terms[i]
                p=rand_path(a, w, t, (used-{w}) | (set(terms)-{t}), rng)
                if p is None: ok=False; break
                legs[i]=p; used |= set(p[1:])
            if ok:
                v=[]
                for i in range(3):
                    p=legs[i]
                    v.append(sum(1 for k in range(len(p)-1)
                                 if frozenset((p[k],p[k+1])) in N) % 2)
                out.add(tuple(v))
        if len(out)==8: break
    return out

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 11)
st=collections.Counter(); misses=[]
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 40):
    a,terms=make_bridge(rng.choice([10,12,14]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    t=tau_edge(E)
    if t<2: continue
    st['bridges with tau>=2']+=1
    best=None
    for trial in range(200):
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
    P={v:y_parities(a,v,terms,N,rng) for v in a if v not in terms}
    fired=False
    for e in E:
        if e in N: continue
        u,v=tuple(e)
        if u in terms or v in terms: continue
        for val in ((0,0,0),(1,1,1)):
            if val in P[u] and val in P[v]: fired=True; break
        if fired: break
    st['PARITY Y-SHIFT FIRES' if fired else '*** still does not fire ***']+=1
    if not fired: misses.append((len(a),len(N),t))
for k in sorted(st): print('  %-40s %s'%(k,st[k]))
if misses: print('  misses (n, |N|, tau):', misses[:8])
