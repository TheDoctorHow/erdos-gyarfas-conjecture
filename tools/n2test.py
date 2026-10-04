"""n2test.py -- try to refute (N2): tau(B) >= 2 ==> (L2**) holds for B.

The threat is the k=3 cactus of 296 with NON-bipartite parts: then no cut edge is
bipartizing (each B-e_i still contains an odd cycle inside a part), so tau >= 2, and the
route structure that killed (L2**) at tau=1 is still present.  Does it still kill it?
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split('"""', 2)[2])

def bipn(a):
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

def relab(E, tag):
    return {frozenset(('%s%s' % (tag, u), '%s%s' % (tag, v))) for u, v in map(tuple, E)}

def evenpaths(a, s, t):
    out = []; p = [s]; on = {s}
    def rec(u):
        for w in a[u]:
            if w == t:
                q = p+[t]
                if (len(q)-1) % 2 == 0: out.append(q)
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out

def es(p): return frozenset(frozenset((p[i], p[i+1])) for i in range(len(p)-1))

def cactus3(ms, rng):
    """three parts in a cycle, each a cubic girth-5 graph minus an edge (NON-bipartite),
    with one further edge subdivided to give the terminal."""
    parts = []
    for k, m in enumerate(ms):
        H = rand_cubic_g5(m, rng)
        if H is None: return None, None
        H = relab(H, 'p%d.' % k)
        e = rng.choice(sorted(H, key=lambda s: sorted(s)))
        u, v = tuple(e)
        parts.append((H - {e}, u, v))
    E = set()
    for (H, u, v) in parts: E |= H
    E |= {frozenset((parts[0][2], parts[1][1])),
          frozenset((parts[1][2], parts[2][1])),
          frozenset((parts[2][2], parts[0][1]))}
    new = set(E); terms = []
    for k in range(3):
        inside = [e for e in new if all(str(x).startswith('p%d.' % k) for x in e)]
        if not inside: return None, None
        e = rng.choice(sorted(inside, key=lambda s: sorted(map(str, s))))
        x, y = tuple(e); t = 'T%d' % k
        new.discard(e); new.add(frozenset((x, t))); new.add(frozenset((y, t)))
        terms.append(t)
    a = adj_of(new)
    if girth(a) < 5 or not two_connected(a) or bipn(a) is not None: return None, None
    return a, sorted(terms)

rng = random.Random(int(sys.argv[1]))
st = collections.Counter()
for _ in range(int(sys.argv[2])):
    a, terms = cactus3([rng.choice([10,12])]*3, rng)
    if a is None: continue
    E = {frozenset((u, w)) for u in a for w in a[u]}
    bip = [e for e in E if bipn(adj_of(E-{e})) is not None]
    st['bridges'] += 1; st['bipartizing_edges_%d' % len(bip)] += 1
    nfail = 0
    for x, y in itertools.combinations(terms, 2):
        ev = evenpaths(a, x, y)
        if len(ev) < 2: nfail += 1; continue
        Es = [es(p) for p in ev]
        nfail += not any(not (Es[i] & Es[j]) for i in range(len(Es))
                         for j in range(i+1, len(Es)))
    st['failing_pairs_%d' % nfail] += 1
    if nfail == 3:
        st['*** (N2) REFUTED ***'] += 1
        print('(N2) COUNTEREXAMPLE n=%d bip=%d' % (len(a), len(bip)))
        print('  edges:', sorted(tuple(sorted(map(str,e))) for e in E))
for k in sorted(st): print('  %-30s %s' % (k, st[k]))
