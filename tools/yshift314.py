"""yshift314.py -- the THREE-TERMINAL mechanism, transplanted to non-bipartite bridges.

Fix a 2-colouring c of V(B) with ALL THREE TERMINALS THE SAME COLOUR; let N be the
non-crossing edges and W = B - N, which is bipartite with colouring c.
If w, w' are ADJACENT in W and each carries a Y - three internally disjoint paths inside W
to the three terminals - then:
   * every leg has mu_m == [c(w) =/= c(b_m)] (mod 2), so every mu_i + mu_j is EVEN, so all
     three concatenations are EVEN connecting paths of B and (SC) pins all three mod 4;
   * subtracting the three relations at w and w' forces every shift e_m EVEN (199.2a);
   * but c(w') =/= c(w), so every e_m is ODD.
Contradiction.  This tests whether the configuration is available in bridges with tau >= 2.
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

def has_Y(a, w, terms):
    """three internally disjoint paths from w to the three terminals, inside a"""
    if w in terms: return False
    # Menger via a simple flow: unit vertex capacities, source w, sink super-terminal
    # brute force for small graphs: try all triples of neighbours as first steps
    nb=sorted(a[w], key=str)
    if len(nb)<3: return False
    def reach(start, banned, tgt):
        if start in banned: return None
        d={start:None}; q=deque([start])
        while q:
            u=q.popleft()
            if u==tgt:
                p=[];z=u
                while z is not None: p.append(z); z=d[z]
                return p
            for x in a[u]:
                if x in banned or x in d: continue
                d[x]=u; q.append(x)
        return None
    for perm in itertools.permutations(terms):
        used={w}; ok=True
        for k,t in enumerate(perm):
            p=None
            for start in nb:
                if start in used: continue
                p=reach(start, used-{start} | {w} , t)
                if p: break
            if not p: ok=False; break
            used |= set(p)
        if ok: return True
    return False

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 5)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 40):
    a,terms=make_bridge(rng.choice([10,12,14]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    t=tau_edge(E)
    if t<2: continue
    st['bridges with tau>=2']+=1
    # colouring: start from a max-cut-ish greedy, then force the terminals to colour 0
    best=None
    for trial in range(300):
        col={v:rng.randint(0,1) for v in a}
        for x in terms: col[x]=0
        improved=True
        while improved:
            improved=False
            for v in a:
                if v in terms: continue
                nc=sum(1 for w in a[v] if col[w]==col[v])
                if nc>len(a[v])-nc: col[v]^=1; improved=True
        N={e for e in E if col[tuple(e)[0]]==col[tuple(e)[1]]}
        if best is None or len(N)<len(best[1]): best=(dict(col),set(N))
    col,N=best
    st['min |N| with terminals monochromatic = %d'%len(N)]+=1
    W=E-N; aW=adj_of(W)
    fired=False
    for u in aW:
        if fired: break
        for v in aW[u]:
            if str(u)>=str(v): continue
            if has_Y(aW,u,terms) and has_Y(aW,v,terms): fired=True; break
    st['Y-SHIFT FIRES' if fired else '*** does not fire ***']+=1
for k in sorted(st): print('  %-42s %s'%(k,st[k]))
