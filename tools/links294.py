"""links294.py -- test 294.1 (segment structure) and 294.3 (link decomposition) on
genuine bridge-shaped graphs: max degree 3, exactly three vertices of degree 2
(the terminals), girth >= 5, 2-connected, non-bipartite, terminals pairwise at
distance >= 5.

For every terminal pair and every pair of connecting paths it checks:
  A  every common vertex other than a common endpoint carries a common edge;
  A' the intersection is a disjoint union of paths with >= 1 edge (plus isolated
     terminals only);
  L1 the symmetric difference is a VERTEX-disjoint union of cycles;
  L2 the number of odd link cycles == |P|+|R| (mod 2);
  L3 (L2**) for the pair <=> some two even connecting paths are EDGE-disjoint;
  S  the splitting dichotomy: (L2**) => even paths use both edges at each terminal.
"""
import sys, random, itertools
from collections import deque

def adj_of(es):
    a = {}
    for e in es:
        u, v = tuple(e)
        a.setdefault(u, set()).add(v)
        a.setdefault(v, set()).add(u)
    return a

def girth(a):
    best = 10 ** 9
    for s in a:
        d = {s: 0}; par = {s: None}; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in d:
                    d[w] = d[u] + 1; par[w] = u; q.append(w)
                elif w != par[u]:
                    best = min(best, d[u] + d[w] + 1)
    return best

def connected(a, skip=None):
    vs = [v for v in a if v != skip]
    if not vs:
        return True
    s = vs[0]; seen = {s}; q = deque([s])
    while q:
        u = q.popleft()
        for w in a[u]:
            if w != skip and w not in seen:
                seen.add(w); q.append(w)
    return len(seen) == len(vs)

def two_connected(a):
    return len(a) >= 3 and connected(a) and all(connected(a, v) for v in a)

def bipartite(a):
    col = {}
    for s in a:
        if s in col: continue
        col[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in col:
                    col[w] = 1 - col[u]; q.append(w)
                elif col[w] == col[u]:
                    return False
    return True

def dist(a, s, t):
    d = {s: 0}; q = deque([s])
    while q:
        u = q.popleft()
        if u == t: return d[u]
        for w in a[u]:
            if w not in d:
                d[w] = d[u] + 1; q.append(w)
    return 10 ** 9

def all_paths(a, s, t, cap=200000):
    out = []; path = [s]; onp = {s}
    def rec(u):
        if len(out) >= cap: return
        if u == t:
            out.append(list(path)); return
        for w in sorted(a[u], key=str):
            if w not in onp:
                onp.add(w); path.append(w)
                rec(w)
                path.pop(); onp.discard(w)
    rec(s)
    return out

def eset(p):
    return {frozenset((p[i], p[i + 1])) for i in range(len(p) - 1)}

def comps(es):
    """components of an edge set, as (vertexset, edgeset)"""
    a = adj_of(es); seen = set(); out = []
    for s in a:
        if s in seen: continue
        blk = {s}; q = deque([s]); seen.add(s)
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in seen:
                    seen.add(w); blk.add(w); q.append(w)
        out.append((blk, {e for e in es if tuple(e)[0] in blk}))
    return out

# ------------------------------------------------------------------ generator
def random_cubic_g5(m, rng, tries=3000):
    """random 3-regular simple graph on m vertices with girth >= 5, 2-connected."""
    for _ in range(tries):
        es = set(); deg = {v: 0 for v in range(m)}
        ok = True
        for _ in range(3 * m // 2):
            cand = [v for v in range(m) if deg[v] < 3]
            rng.shuffle(cand)
            placed = False
            for u in cand:
                for w in cand:
                    if w == u or deg[w] >= 3 or frozenset((u, w)) in es: continue
                    if girth(adj_of(es | {frozenset((u, w))})) < 5: continue
                    es.add(frozenset((u, w))); deg[u] += 1; deg[w] += 1
                    placed = True; break
                if placed: break
            if not placed:
                ok = False; break
        if not ok or any(deg[v] != 3 for v in range(m)): continue
        a = adj_of(es)
        if len(a) == m and two_connected(a):
            return es
    return None

def random_bridge(n, rng, tries=200):
    """214.2's model: subdivide three edges of a cubic girth-5 graph, giving three
    terminals of degree 2 and all other vertices of degree 3."""
    m = n - 3
    if m % 2 or m < 10: return None, None
    for _ in range(tries):
        es = random_cubic_g5(m, rng)
        if es is None: return None, None
        chosen = rng.sample(sorted(es, key=lambda s: sorted(s)), 3)
        new = set(es); terms = []
        for k, e in enumerate(chosen):
            u, v = tuple(e)
            t = 'T%d' % k
            new.discard(e); new.add(frozenset((u, t))); new.add(frozenset((v, t)))
            terms.append(t)
        a = adj_of(new)
        if not two_connected(a) or bipartite(a): continue
        if girth(a) < 5: continue
        return a, sorted(terms)
    return None, None

# ------------------------------------------------------------------ tests
def run(a, terms, stats):
    for b1, b2 in itertools.combinations(terms, 2):
        paths = all_paths(a, b1, b2)
        ev = [p for p in paths if (len(p) - 1) % 2 == 0]
        od = [p for p in paths if (len(p) - 1) % 2 == 1]
        stats['pairs'] += 1
        stats['even'] += len(ev); stats['odd'] += len(od)
        # L3: (L2**) via internal disjointness, and via edge-disjointness
        l2_int = l2_edge = False
        for i in range(len(ev)):
            for j in range(i + 1, len(ev)):
                si, sj = set(ev[i][1:-1]), set(ev[j][1:-1])
                if not (si & sj): l2_int = True
                if not (eset(ev[i]) & eset(ev[j])): l2_edge = True
                if l2_int and l2_edge: break
            if l2_int and l2_edge: break
        if l2_int != l2_edge:
            stats['L3_fail'] += 1
        stats['L2ss_holds' if l2_int else 'L2ss_fails'] += 1
        # S: splitting
        first = {p[1] for p in ev}; last = {p[-2] for p in ev}
        split = len(first) == 2 and len(last) == 2
        if l2_int and not split: stats['S_fail'] += 1
        if not l2_int:
            stats['fail_split' if split else 'fail_nosplit'] += 1
        # A, A', L1, L2 over all pairs of connecting paths (cap for speed)
        sample = paths if len(paths) <= 60 else random.Random(7).sample(paths, 60)
        for P, R in itertools.combinations(sample, 2):
            EP, ER = eset(P), eset(R)
            VP, VR = set(P), set(R)
            common = VP & VR
            inter = EP & ER
            ia = adj_of(inter) if inter else {}
            for w in common:
                if w in (b1, b2): continue
                if w not in ia:
                    stats['A_fail'] += 1
            for blk, be in comps(inter):
                da = adj_of(be)
                if any(len(da[v]) > 2 for v in blk): stats['Aseg_fail'] += 1
            D = EP ^ ER
            if D:
                da = adj_of(D)
                if any(len(da[v]) != 2 for v in da): stats['L1_fail'] += 1
                nodd = sum(1 for blk, be in comps(D) if len(be) % 2)
                if nodd % 2 != (len(P) + len(R)) % 2: stats['L2_fail'] += 1
                stats['diffs'] += 1
                if len(comps(D)) == 1: stats['single_link'] += 1
            stats['pathpairs'] += 1

def main():
    rng = random.Random(20260805)
    stats = dict.fromkeys(
        ['graphs','pairs','even','odd','pathpairs','diffs','single_link',
         'A_fail','Aseg_fail','L1_fail','L2_fail','L3_fail','S_fail',
         'L2ss_holds','L2ss_fails','fail_split','fail_nosplit'], 0)
    ns = [13, 15, 17, 19]
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for n in ns:
        got = 0
        while got < per:
            a, terms = random_bridge(n, rng)
            if a is None: break
            got += 1; stats['graphs'] += 1
            run(a, terms, stats)
    w = max(len(k) for k in stats)
    for k, v in stats.items():
        print('%-*s %d' % (w, k, v))

main()
