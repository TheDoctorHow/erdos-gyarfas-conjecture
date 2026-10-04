"""EXHAUSTIVE TEST of the BINARY EAR CONJECTURE (section 170.6).

A graph G satisfies the conjecture if there is a cycle C and ears P_1..P_k with
 - interiors pairwise disjoint and disjoint from C,
 - chosen arcs A_1..A_k pairwise EDGE-disjoint sub-arcs of C,
 - v_2(delta_j) = j-1 for j=1..k where delta_j = |P_j| - |A_j|,
 - |C| + sum over {j : delta_j>0} of delta_j  <  2^{k+1}.
Then (170.2) G has a cycle of length exactly 2^k.
"""
import itertools, sys

def v2(d):
    if d == 0: return None          # valuation infinite - unusable
    k = 0
    d = abs(d)
    while d % 2 == 0: d //= 2; k += 1
    return k

def from_edges(n, E):
    adj = [[] for _ in range(n)]
    for a, b in E: adj[a].append(b); adj[b].append(a)
    return adj

def all_cycles(adj):
    n = len(adj); out = []
    for s in range(n):
        stack = [(s, [s], 1 << s)]
        while stack:
            u, path, mask = stack.pop()
            for v in adj[u]:
                if v == s and len(path) >= 3:
                    if path[1] < path[-1]: out.append(path[:])
                elif v > s and not (mask >> v) & 1:
                    stack.append((v, path + [v], mask | (1 << v)))
    return out

def ears_of(adj, cyc):
    """all ears: (u, v, length, frozenset(interior))"""
    n = len(adj); VC = set(cyc); L = len(cyc)
    cyc_edges = set()
    for i in range(L):
        a, b = cyc[i], cyc[(i+1) % L]; cyc_edges.add((min(a,b), max(a,b)))
    res = []
    # chords
    for u in cyc:
        for v in adj[u]:
            if v in VC and u < v and (u, v) not in cyc_edges:
                res.append((u, v, 1, frozenset()))
    # longer ears through V \ VC
    outside = [x for x in range(n) if x not in VC]
    for u in cyc:
        starts = [w for w in adj[u] if w not in VC]
        for w in starts:
            stack = [(w, [w], {w})]
            while stack:
                x, path, seen = stack.pop()
                for y in adj[x]:
                    if y in VC:
                        if y != u and u < y:
                            res.append((u, y, len(path)+1, frozenset(path)))
                    elif y not in seen:
                        stack.append((y, path + [y], seen | {y}))
    # dedupe
    return list({(u,v,l,i) for (u,v,l,i) in res})

def arc_edges(cyc, u, v, which):
    """edge set of one of the two arcs between u and v; which in {0,1}"""
    L = len(cyc); pos = {w: i for i, w in enumerate(cyc)}
    p, q = pos[u], pos[v]
    if which == 0:
        idx = list(range(p, q)) if p < q else list(range(p, L)) + list(range(0, q))
    else:
        idx = list(range(q, p)) if q < p else list(range(q, L)) + list(range(0, p))
    es = frozenset((min(cyc[i], cyc[(i+1) % L]), max(cyc[i], cyc[(i+1) % L])) for i in idx)
    return es

def check(adj, name, kmax=6, verbose=False):
    best = None
    for cyc in all_cycles(adj):
        L = len(cyc)
        opts = []          # (delta, v2, arcEdgeSet, interior)
        for (u, v, ln, inter) in ears_of(adj, cyc):
            for which in (0, 1):
                ae = arc_edges(cyc, u, v, which)
                d = ln - len(ae)
                val = v2(d)
                if val is None: continue
                opts.append((d, val, ae, inter))
        # group by valuation
        byval = {}
        for o in opts: byval.setdefault(o[1], []).append(o)
        for k in range(2, kmax+1):
            if any(j not in byval for j in range(k)): continue
            for combo in itertools.product(*[byval[j] for j in range(k)]):
                ok = True
                for a, b in itertools.combinations(combo, 2):
                    if a[2] & b[2]: ok = False; break
                    if a[3] & b[3]: ok = False; break
                if not ok: continue
                pos_sum = sum(d for (d, _, _, _) in combo if d > 0)
                if L + pos_sum < 2**(k+1):
                    if best is None or k > best[0]:
                        best = (k, L, [d for (d,_,_,_) in combo])
    if best:
        print("%-16s SATISFIED: k=%d, |C|=%d, deltas=%s  -> cycle of length %d"
              % (name, best[0], best[1], best[2], 2**best[0]))
        return True
    else:
        print("%-16s *** NO BINARY EAR BASIS FOUND ***" % name)
        return False

PET = from_edges(10, [(0,1),(1,2),(2,3),(3,4),(4,0),
                      (0,5),(1,6),(2,7),(3,8),(4,9),
                      (5,7),(7,9),(9,6),(6,8),(8,5)])
K33 = from_edges(6, [(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)])
K4  = from_edges(4, [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)])
Q3  = from_edges(8, [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
                     (0,4),(1,5),(2,6),(3,7)])
def lcf(n, sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return from_edges(n, sorted({tuple(sorted(e)) for e in E}))
HEAWOOD = lcf(14, [5,-5])
MK      = lcf(16, [5,-5])

print("BINARY EAR CONJECTURE -- exhaustive check\n")
for nm, g in [("K4", K4), ("K3,3", K33), ("Q3 (cube)", Q3), ("PETERSEN", PET),
              ("Heawood", HEAWOOD), ("Mobius-Kantor", MK)]:
    check(adj=g, name=nm)
    sys.stdout.flush()
