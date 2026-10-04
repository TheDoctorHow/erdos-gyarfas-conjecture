"""ecore2.py -- try to BREAK (E-CORE) at larger n.

(E-CORE)  a terminal pair FAILS (L2**)  <=>  some EDGE lies on EVERY even connecting
          path between them.

Population as in 214.2: configuration-model 3-regular graph on m vertices, rejected
unless simple, girth >= 5 and 2-connected; then three edges subdivided to make the
terminals.  Reports every failing pair and whether it has a common edge.
"""
import sys, random, itertools, collections
from collections import deque

def adj_of(es):
    a = {}
    for e in es:
        u, v = tuple(e)
        a.setdefault(u, set()).add(v); a.setdefault(v, set()).add(u)
    return a

def girth(a):
    best = 10 ** 9
    for s in a:
        d = {s: 0}; par = {s: None}; q = deque([s])
        while q:
            u = q.popleft()
            if d[u] * 2 >= best: continue
            for w in a[u]:
                if w not in d:
                    d[w] = d[u] + 1; par[w] = u; q.append(w)
                elif w != par[u]:
                    best = min(best, d[u] + d[w] + 1)
    return best

def connected(a, skip=None):
    vs = [v for v in a if v != skip]
    if not vs: return True
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
                if w not in col: col[w] = 1 - col[u]; q.append(w)
                elif col[w] == col[u]: return False
    return True

def rand_cubic_g5(m, rng, tries=20000):
    for _ in range(tries):
        stubs = [v for v in range(m) for _ in range(3)]
        rng.shuffle(stubs)
        es = set(); ok = True
        for i in range(0, len(stubs), 2):
            u, v = stubs[i], stubs[i + 1]
            if u == v or frozenset((u, v)) in es: ok = False; break
            es.add(frozenset((u, v)))
        if not ok or len(es) != 3 * m // 2: continue
        a = adj_of(es)
        if len(a) != m or girth(a) < 5 or not two_connected(a): continue
        return es
    return None

def make_bridge(m, rng):
    es = rand_cubic_g5(m, rng)
    if es is None: return None, None
    chosen = rng.sample(sorted(es, key=lambda s: sorted(s)), 3)
    new = set(es); terms = []
    for k, e in enumerate(chosen):
        u, v = tuple(e); t = 'T%d' % k
        new.discard(e); new.add(frozenset((u, t))); new.add(frozenset((v, t)))
        terms.append(t)
    a = adj_of(new)
    if not two_connected(a) or bipartite(a) or girth(a) < 5: return None, None
    return a, sorted(terms)

def even_paths(a, s, t):
    """all s-t paths, returned as (vertexlist, edgeframeset) for EVEN ones only"""
    out = []; path = [s]; onp = {s}
    def rec(u):
        for w in a[u]:
            if w == t:
                if len(path) % 2 == 1:          # length = len(path), even iff len(path)%2==0...
                    pass
                p = path + [t]
                if (len(p) - 1) % 2 == 0:
                    out.append(p)
            elif w not in onp:
                onp.add(w); path.append(w); rec(w); path.pop(); onp.discard(w)
    rec(s)
    return out

def eset(p):
    return frozenset(frozenset((p[i], p[i + 1])) for i in range(len(p) - 1))

def main():
    seed = int(sys.argv[1]); mlist = [int(x) for x in sys.argv[2].split(',')]
    per = int(sys.argv[3])
    rng = random.Random(seed)
    st = collections.Counter()
    for m in mlist:
        got = 0; attempts = 0
        while got < per and attempts < per * 200:
            attempts += 1
            a, terms = make_bridge(m, rng)
            if a is None: continue
            got += 1; st['bridges'] += 1
            for b1, b2 in itertools.combinations(terms, 2):
                ev = even_paths(a, b1, b2)
                if len(ev) < 2:
                    st['degenerate'] += 1; continue
                E = [eset(p) for p in ev]
                st['pairs'] += 1
                st['maxpaths'] = max(st['maxpaths'], len(ev))
                # (E-CORE) easy direction: a common edge forces failure.  So compute the
                # core first; only when it is EMPTY need we search for a disjoint pair,
                # and (E-CORE) predicts we will always find one.
                core = set(E[0])
                for x in E[1:]: core &= x
                if core:
                    st['FAIL'] += 1
                    st['fail_with_common_edge'] += 1
                    st['core_%d' % len(core)] += 1
                    continue
                idx = list(range(len(E))); rng.shuffle(idx)
                found = False
                for ii in range(len(idx)):
                    Ei = E[idx[ii]]
                    for jj in range(ii + 1, len(idx)):
                        if not (Ei & E[idx[jj]]): found = True; break
                    if found: break
                if found:
                    st['HOLDS'] += 1
                else:
                    st['FAIL'] += 1
                    st['*** E-CORE COUNTEREXAMPLE ***'] += 1
                    print('COUNTEREXAMPLE  n=%d  terminals %s %s  #even=%d'
                          % (len(a), b1, b2, len(ev)))
                    print('  edges:', sorted(tuple(sorted(e)) for e in
                                             {x for p in ev for x in eset(p)}))
                    print('  graph:', sorted(tuple(sorted(map(str, e))) for e in
                                             {frozenset((u, w)) for u in a for w in a[u]}))
    print('m=%s  seed=%d' % (mlist, seed))
    for k in sorted(st): print('  %-32s %s' % (k, st[k]))

main()
