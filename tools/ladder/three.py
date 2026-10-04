"""three.py -- can ALL THREE terminal pairs of a bridge fail (L2**)?

(L2**) only asks for SOME pair, so this is the statement that really decides it.
250.4 and 251.1 measured 0 of 280 configurations with all three pairs failing.

Predicted configuration (from the parity identity of 295.1): three parts D_1, D_2, D_3
each holding one terminal, joined in a CYCLE by three single edges e_1, e_2, e_3.  Then
every b_i-b_j path uses either one cut edge or the other two, so for every pair the even
paths all share an edge and no two are edge-disjoint - all three pairs fail.

Such a B has three 2-edge-cuts, so by 295.8 it exhibits nontrivial 3-edge-cuts of G:
it can only occur if the minimal counterexample is NOT cyclically 4-edge-connected.
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

def paths(a, s, t):
    out = []; p = [s]; on = {s}
    def rec(u):
        for w in a[u]:
            if w == t: out.append(p + [t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out

def es(p): return frozenset(frozenset((p[i], p[i+1])) for i in range(len(p)-1))

def relabel(E, tag):
    return {frozenset(('%s%s' % (tag, u), '%s%s' % (tag, v))) for u, v in map(tuple, E)}

def build(ms, rng):
    """three cubic girth-5 parts, one edge removed from each, joined in a cycle by three
    edges; then one edge of each part subdivided to give one terminal per part."""
    parts = []
    for k, m in enumerate(ms):
        H = rand_cubic_g5(m, rng)
        if H is None: return None, None
        H = relabel(H, 'p%d' % k)
        e = rng.choice(sorted(H, key=lambda s: sorted(s)))
        u, v = tuple(e)
        parts.append((H - {e}, u, v))
    E = set()
    for (H, u, v) in parts: E |= H
    # join in a cycle: part0.v - part1.u , part1.v - part2.u , part2.v - part0.u
    E |= {frozenset((parts[0][2], parts[1][1])),
          frozenset((parts[1][2], parts[2][1])),
          frozenset((parts[2][2], parts[0][1]))}
    a = adj_of(E)
    if girth(a) < 5 or not two_connected(a): return None, None
    terms = []
    new = set(E)
    for k, (H, u, v) in enumerate(parts):
        inside = [e for e in new if all(str(x).startswith('p%d' % k) for x in e)]
        if not inside: return None, None
        e = rng.choice(sorted(inside, key=lambda s: sorted(map(str, s))))
        x, y = tuple(e); t = 'T%d' % k
        new.discard(e); new.add(frozenset((x, t))); new.add(frozenset((y, t)))
        terms.append(t)
    a = adj_of(new)
    if girth(a) < 5 or not two_connected(a) or bipartition(a) is not None:
        return None, None
    def d(s0, t0):
        D = {s0: 0}; q = deque([s0])
        while q:
            u0 = q.popleft()
            if u0 == t0: return D[u0]
            for w in a[u0]:
                if w not in D: D[w] = D[u0] + 1; q.append(w)
        return 10**9
    dist = min(d(x, y) for x, y in itertools.combinations(terms, 2))
    return a, (sorted(terms), dist)

def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
    tries = int(sys.argv[2]) if len(sys.argv) > 2 else 25
    st = collections.Counter(); per = collections.Counter(); shown = 0
    for _ in range(tries):
        a, info = build([rng.choice([10, 12])] * 3, rng)
        if a is None: continue
        terms, dist = info
        st['bridges'] += 1; st['min_terminal_distance_%d' % dist] += 1
        E = {frozenset((u, w)) for u in a for w in a[u]}
        nfail = 0; cores = []
        for b1, b2 in itertools.combinations(terms, 2):
            ev = [p for p in paths(a, b1, b2) if (len(p) - 1) % 2 == 0]
            if len(ev) < 2: st['degenerate'] += 1; continue
            Es = [es(p) for p in ev]
            core = set(Es[0])
            for x in Es[1:]: core &= x
            fail = bool(core) or not any(not (Es[i] & Es[j])
                        for i in range(len(Es)) for j in range(i + 1, len(Es)))
            st['pairs'] += 1; st['FAIL' if fail else 'HOLDS'] += 1
            if fail:
                nfail += 1
                cores.append((b1, b2, sorted(tuple(sorted(map(str, e))) for e in core)))
        per[nfail] += 1
        if nfail == 3 and shown < 2:
            shown += 1
            st['*** ALL THREE PAIRS FAIL ***'] += 1
            print('ALL-THREE-FAIL  n=%d  min terminal distance %d' % (len(a), dist))
            for c in cores: print('   pair %s,%s  common edges %s' % c)
            print('   edges:', sorted(tuple(sorted(map(str, e))) for e in E))
    for k in sorted(st): print('  %-34s %s' % (k, st[k]))
    print('  failing pairs per bridge:', dict(sorted(per.items())))

main()
