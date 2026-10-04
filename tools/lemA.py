"""lemA.py -- LEMMA A: in a 2-connected graph, EVERY cycle is the symmetric difference
of two s-t paths (for any s =/= t).  This is strictly stronger than (SPAN)."""
import itertools, random, collections
from collections import deque

def adj(E):
    a = {}
    for e in E:
        u, v = tuple(e); a.setdefault(u, set()).add(v); a.setdefault(v, set()).add(u)
    return a
def conn(a, skip=None):
    vs = [v for v in a if v != skip]
    if not vs: return True
    s = vs[0]; seen = {s}; q = deque([s])
    while q:
        u = q.popleft()
        for w in a[u]:
            if w != skip and w not in seen: seen.add(w); q.append(w)
    return len(seen) == len(vs)
def two_conn(a, n): return len(a) == n and conn(a) and all(conn(a, v) for v in a)
def paths(a, s, t):
    out = []; p = [s]; on = {s}
    def rec(u):
        for w in a[u]:
            if w == t: out.append(list(p)+[t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out
def eset(p): return frozenset(frozenset((p[i], p[i+1])) for i in range(len(p)-1))
def cycles(E):
    """all cycles as edge sets"""
    a = adj(E); out = set(); vs = sorted(a, key=str); idx = {v: i for i, v in enumerate(vs)}
    for i, s in enumerate(vs):
        p = [s]; on = {s}
        def rec(u):
            for w in a[u]:
                if idx[w] < i: continue
                if w == s and len(p) >= 3: out.add(eset(p+[s]))
                elif w not in on:
                    on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
        rec(s)
    return out

def check(E, n, st, label):
    a = adj(E)
    if not two_conn(a, n): return
    Cy = cycles(E)
    for s, t in itertools.combinations(sorted(a, key=str), 2):
        P = [eset(p) for p in paths(a, s, t)]
        D = {p1 ^ p2 for p1 in P for p2 in P}
        for Z in Cy:
            st['cycle,pair tested'] += 1
            if Z in D: st['is a difference'] += 1
            else:
                st['*** LEMMA A FAILS ***'] += 1
                if st['*** LEMMA A FAILS ***'] <= 2:
                    print('FAIL %s  s=%s t=%s  Z=%s' % (label, s, t,
                          sorted(tuple(sorted(map(str, e))) for e in Z)))

st = collections.Counter()
for n in range(3, 7):
    allE = [frozenset(p) for p in itertools.combinations(range(n), 2)]
    for k in range(n, len(allE)+1):
        for S in itertools.combinations(allE, k):
            check(set(S), n, st, 'exh n=%d' % n)
    print('  ...n=%d done: %s' % (n, dict(st)))
rng = random.Random(8)
for _ in range(150):
    n = rng.choice([7, 8, 9, 10])
    allE = [frozenset(p) for p in itertools.combinations(range(n), 2)]
    check(set(rng.sample(allE, rng.randint(n, min(len(allE), 2*n)))), n, st, 'rand')
for k in sorted(st): print('  %-24s %s' % (k, st[k]))
