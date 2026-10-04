"""twocut.py -- can a bridge have TWO bipartizing edges, and can all THREE terminal
pairs then fail?

Two distinct edges e,f with B-e and B-f both bipartite necessarily form a 2-EDGE-CUT of
B (proved: otherwise B-e-f is connected and the two colourings coincide, forcing f both
to cross and not to cross).  Random cubic graphs are 3-edge-connected, so the sweep so
far could never produce one.  Here they are built deliberately: take two cubic girth-5
graphs, delete an edge from each, and join the four resulting degree-2 vertices by two
edges; then subdivide three edges to make the terminals.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split('"""', 2)[2])

def bipartition(a):
    col = {}
    for s in a:
        if s in col: continue
        col[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in col: col[w] = 1 - col[u]; q.append(w)
                elif col[w] == col[u]: return None
    return col

def even_paths_all(a, s, t):
    out = []; path = [s]; onp = {s}
    def rec(u):
        for w in a[u]:
            if w == t:
                p = path + [t]
                if (len(p) - 1) % 2 == 0: out.append(p)
            elif w not in onp:
                onp.add(w); path.append(w); rec(w); path.pop(); onp.discard(w)
    rec(s); return out

def relabel(es, tag):
    return {frozenset(('%s%s' % (tag, u), '%s%s' % (tag, v))) for u, v in
            (tuple(e) for e in es)}

def build(m1, m2, rng):
    H1 = rand_cubic_g5(m1, rng); H2 = rand_cubic_g5(m2, rng)
    if H1 is None or H2 is None: return None, None
    H1 = relabel(H1, 'a'); H2 = relabel(H2, 'b')
    e1 = rng.choice(sorted(H1, key=lambda s: sorted(s)))
    e2 = rng.choice(sorted(H2, key=lambda s: sorted(s)))
    p1, q1 = tuple(e1); p2, q2 = tuple(e2)
    E = (H1 - {e1}) | (H2 - {e2}) | {frozenset((p1, p2)), frozenset((q1, q2))}
    a = adj_of(E)
    if girth(a) < 5 or not two_connected(a): return None, None
    # subdivide three edges: one inside H1's side, two inside H2's side, so the
    # 2-edge-cut has one terminal on one side and two on the other
    c1 = [e for e in E if all(str(x).startswith('a') for x in e)]
    c2 = [e for e in E if all(str(x).startswith('b') for x in e)]
    if len(c1) < 1 or len(c2) < 2: return None, None
    pick = [rng.choice(c1)] + rng.sample(c2, 2)
    new = set(E); terms = []
    for k, e in enumerate(pick):
        u, v = tuple(e); t = 'T%d' % k
        new.discard(e); new.add(frozenset((u, t))); new.add(frozenset((v, t)))
        terms.append(t)
    a = adj_of(new)
    if girth(a) < 5 or not two_connected(a) or bipartition(a) is not None:
        return None, None
    return a, sorted(terms)

def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
    tries = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    st = collections.Counter(); perbridge = collections.Counter()
    for _ in range(tries):
        a, terms = build(rng.choice([10, 12]), rng.choice([10, 12]), rng)
        if a is None: continue
        st['bridges'] += 1
        E = {frozenset((u, w)) for u in a for w in a[u]}
        bipedges = [e for e in E if bipartition(adj_of(E - {e})) is not None]
        st['bipartizing_edges_%d' % len(bipedges)] += 1
        nfail = 0; witnessed = 0
        for b1, b2 in itertools.combinations(terms, 2):
            ev = even_paths_all(a, b1, b2)
            if len(ev) < 2: st['degenerate'] += 1; continue
            Es = [frozenset(frozenset((p[i], p[i+1])) for i in range(len(p)-1))
                  for p in ev]
            core = set(Es[0])
            for x in Es[1:]: core &= x
            fail = bool(core) or not any(not (Es[i] & Es[j])
                        for i in range(len(Es)) for j in range(i+1, len(Es)))
            st['pairs'] += 1
            if not fail: st['HOLDS'] += 1; continue
            nfail += 1; st['FAIL'] += 1
            w = [e for e in bipedges
                 if bipartition(adj_of(E - {e}))[b1] != bipartition(adj_of(E - {e}))[b2]]
            if w: witnessed += 1; st['fail_with_BIP_witness'] += 1
            else: st['BIP(sharpened) has NO witness'] += 1
            if core:
                st['ECORE_ok'] += 1; st['core_size_%d' % len(core)] += 1
            else:
                st['*** E-CORE COUNTEREXAMPLE ***'] += 1
                print('E-CORE COUNTEREXAMPLE n=%d pair=%s,%s #even=%d'
                      % (len(a), b1, b2, len(ev)))
                print('   edges:', sorted(tuple(sorted(map(str, x))) for x in E))
        perbridge[nfail] += 1
        if nfail == 3:
            st['ALL THREE PAIRS FAIL'] += 1
            print('ALL-THREE-FAIL bridge, n=%d, bipartizing edges=%s'
                  % (len(a), sorted(tuple(sorted(map(str, e))) for e in bipedges)))
            print('  edges:', sorted(tuple(sorted(map(str, e))) for e in E))
    for k in sorted(st): print('  %-32s %s' % (k, st[k]))
    print('  failing pairs per bridge:', dict(sorted(perbridge.items())))

main()
