"""l2fail.py -- an explicit bridge-shaped graph for which (L2**) FAILS FOR ALL THREE
TERMINAL PAIRS.

Construction (predicted by the parity identity 295.1):
  D = Heawood graph minus one vertex.  Heawood is cubic, bipartite, girth 6, with parts
  of size 7; deleting a vertex of one part leaves 13 vertices, 18 edges, girth >= 6,
  2-connected, and EXACTLY THREE vertices of degree 2 -- the neighbours of the deleted
  vertex -- all in the SAME class.  Use one as a terminal and the other two as the ends
  of the cut edges.
  B = three copies of D joined in a CYCLE by three edges.

Then |Bip(B)| = 3, and by 295.1 every terminal pair's even connecting paths all run
through a common cut edge, so no two are edge-disjoint: all three pairs fail.
Everything below is checked directly; the failure of each pair is certified WITHOUT
enumerating paths, by the parity identity plus the route structure.
"""
import itertools
from collections import deque

def heawood():
    E = set()
    for i in range(14): E.add(frozenset((i, (i + 1) % 14)))
    for i in range(0, 14, 2): E.add(frozenset((i, (i + 5) % 14)))
    return E

def adj(E):
    a = {}
    for e in E:
        u, v = tuple(e); a.setdefault(u, set()).add(v); a.setdefault(v, set()).add(u)
    return a

def girth(a):
    best = 10**9
    for s in a:
        d = {s: 0}; par = {s: None}; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in d: d[w] = d[u]+1; par[w] = u; q.append(w)
                elif w != par[u]: best = min(best, d[u]+d[w]+1)
    return best

def conn(a, skipv=None, skipE=frozenset()):
    vs = [v for v in a if v != skipv]
    if not vs: return True
    s = vs[0]; seen = {s}; q = deque([s])
    while q:
        u = q.popleft()
        for w in a[u]:
            if w == skipv or w in seen: continue
            if frozenset((u, w)) in skipE: continue
            seen.add(w); q.append(w)
    return len(seen) == len(vs)

def two_conn(a): return all(conn(a, skipv=v) for v in a)

def bip(a):
    col = {}
    for s in a:
        if s in col: continue
        col[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in col: col[w] = 1-col[u]; q.append(w)
                elif col[w] == col[u]: return None
    return col

def dist(a, s, t):
    d = {s: 0}; q = deque([s])
    while q:
        u = q.popleft()
        if u == t: return d[u]
        for w in a[u]:
            if w not in d: d[w] = d[u]+1; q.append(w)
    return 10**9

# ---- build one part: Heawood minus vertex 13 -----------------------------------
H = heawood()
DEL = 13
NB = sorted({x for e in H if DEL in e for x in e} - {DEL})
D_edges = {e for e in H if DEL not in e}
print('part D = Heawood - v%d :  %d vertices, %d edges, degree-2 vertices %s'
      % (DEL, len({x for e in D_edges for x in e}), len(D_edges), NB))
aD = adj(D_edges)
print('   girth %d, 2-connected %s, bipartite %s, deg-2 all in one class %s'
      % (girth(aD), two_conn(aD), bip(aD) is not None,
         len({bip(aD)[x] for x in NB}) == 1))

T_LOCAL, U_LOCAL, V_LOCAL = NB          # terminal, and the two cut-endpoints

# ---- glue three copies in a cycle ---------------------------------------------
E = set(); T = []; U = []; V = []
for k in range(3):
    for e in D_edges:
        u, v = tuple(e); E.add(frozenset(('%d.%d' % (k, u), '%d.%d' % (k, v))))
    T.append('%d.%d' % (k, T_LOCAL)); U.append('%d.%d' % (k, U_LOCAL))
    V.append('%d.%d' % (k, V_LOCAL))
CUT = [frozenset((V[k], U[(k+1) % 3])) for k in range(3)]
E |= set(CUT)
a = adj(E)

print()
print('B: n = %d, m = %d, max degree %d' % (len(a), len(E), max(len(a[v]) for v in a)))
print('   degree-2 vertices :', sorted(v for v in a if len(a[v]) == 2))
print('   girth %d, 2-connected %s, bipartite %s'
      % (girth(a), two_conn(a), bip(a) is not None))
print('   terminal distances:',
      {(x, y): dist(a, x, y) for x, y in itertools.combinations(T, 2)})
print('   cut edges         :', [tuple(sorted(e)) for e in CUT])

# ---- the bipartizing edges ------------------------------------------------------
Bip = [e for e in E if bip(adj(E - {e})) is not None]
print('   bipartizing edges :', sorted(tuple(sorted(e)) for e in Bip))

# ---- certify the failure of every pair, via 295.1 ------------------------------
print()
print('CERTIFICATE (parity identity 295.1, no path enumeration):')
e1 = Bip[0]
c = bip(adj(E - {e1}))                      # unique 2-colouring of B - e1
N = {e1}                                    # its non-crossing set
for x, y in itertools.combinations(T, 2):
    delta = int(c[x] != c[y])
    # route structure: with the three cut edges removed the parts are separate, so every
    # x-y path uses either the single cut edge between their parts, or the other two.
    kx = T.index(x); ky = T.index(y)
    direct = [e for e in CUT if e == frozenset((V[min(kx, ky)], U[max(kx, ky)]))
              or e == frozenset((V[max(kx, ky)], U[min(kx, ky)]))]
    other = [e for e in CUT if e not in direct]
    uses_e1_direct = int(e1 in direct)
    # a path is EVEN iff |P n N| == delta ; |P n N| is 1 if it uses e1, else 0
    even_route = 'direct' if uses_e1_direct == delta else 'long'
    shared = direct if even_route == 'direct' else other
    print('  pair %s,%s : delta=%d, even paths are the %-6s route, which uses %s'
          % (x, y, delta, even_route,
             ' and '.join(str(tuple(sorted(s))) for s in shared)))
    print('        => every even connecting path contains %s, so no two are edge-'
          'disjoint: THE PAIR FAILS' % str(tuple(sorted(shared[0]))))

print()
print('Each part is separated from the rest of G by its two cut edges plus one')
print('attachment edge: three NONTRIVIAL 3-EDGE-CUTS, each with cycles on both sides.')
print('So this B can only be a bridge if G is NOT cyclically 4-edge-connected.')
