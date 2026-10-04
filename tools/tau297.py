"""tau297.py -- does the frustration index tau stratify the failure of (L2**)?
   Prediction: tau >= 2  ==>  (L2**) holds for B (some pair has two disjoint even paths).
   Population: both the 3-edge-connected bridges (subdivide a random cubic graph) and the
   2-edge-cut bridges (two cubic parts joined by two edges), so the sweep can see the
   configurations that broke (E-CORE)."""
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

def tau(E, cap=2):
    """0 if bipartite, 1 if some single edge bipartizes, else 2 (meaning >= 2)."""
    if bipn(adj_of(set(E))) is not None: return 0
    for e in E:
        if bipn(adj_of(set(E)-{e})) is not None: return 1
    return 2

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

def relab(E, tag):
    return {frozenset(('%s%s' % (tag, u), '%s%s' % (tag, v))) for u, v in map(tuple, E)}

def build2cut(m1, m2, rng):
    H1 = rand_cubic_g5(m1, rng); H2 = rand_cubic_g5(m2, rng)
    if H1 is None or H2 is None: return None, None
    H1 = relab(H1, 'a'); H2 = relab(H2, 'b')
    e1 = rng.choice(sorted(H1, key=lambda s: sorted(s)))
    e2 = rng.choice(sorted(H2, key=lambda s: sorted(s)))
    p1, q1 = tuple(e1); p2, q2 = tuple(e2)
    E = (H1-{e1}) | (H2-{e2}) | {frozenset((p1, p2)), frozenset((q1, q2))}
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
    if girth(a) < 5 or not two_connected(a) or bipn(a) is not None: return None, None
    return a, sorted(terms)

rng = random.Random(int(sys.argv[1]))
st = collections.Counter()
for trial in range(int(sys.argv[2])):
    for maker in (lambda: make_bridge(rng.choice([10,12,14,16]), rng),
                  lambda: build2cut(rng.choice([10,12]), rng.choice([10,12]), rng)):
        a, terms = maker()
        if a is None: continue
        E = {frozenset((u, w)) for u in a for w in a[u]}
        t = tau(E); st['bridges'] += 1; st['tau=%d' % t] += 1
        nfail = 0; npair = 0
        for x, y in itertools.combinations(terms, 2):
            ev = evenpaths(a, x, y)
            if len(ev) < 2: nfail += 1; npair += 1; continue
            Es = [es(p) for p in ev]
            f = not any(not (Es[i] & Es[j]) for i in range(len(Es))
                        for j in range(i+1, len(Es)))
            nfail += f; npair += 1
        st['tau=%d,fail=%d' % (t, nfail)] += 1
        if nfail == npair:
            st['*** (L2**) FAILS, tau=%d ***' % t] += 1
for k in sorted(st): print('  %-28s %s' % (k, st[k]))
