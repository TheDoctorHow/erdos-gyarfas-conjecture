"""Ear / bridge analysis of a longest cycle, to test the Ear Subset-Sum Conjecture.

Setup.  G cubic, C a longest cycle, c = |C|.  An EAR is a path P with both endpoints on
C and all interior vertices off C (a chord is an ear of length 1).  If P joins u,v and
the two arcs of C between them have lengths a and c-a, then replacing an arc by P gives
cycles of lengths (c-a)+|P| and a+|P|, so the two DECREMENTS are

        delta = a - |P|      and      delta' = (c-a) - |P|.

Because C is longest, both are >= 0, i.e. |P| <= min(a, c-a).

Two ears may be applied SIMULTANEOUSLY iff they are internally disjoint and the arcs
they replace are disjoint; then the decrements add.

CONJECTURE UNDER TEST: some longest cycle admits an independent family of ears with
every delta <= 2 and sum of deltas >= c/2.
"""
import sys, itertools, collections
sys.path.insert(0, '.')
from ref import from_graph6


def adjacency(n, edges):
    a = collections.defaultdict(list)
    for u, v in edges:
        a[u].append(v); a[v].append(u)
    return a


def all_cycles(n, adj):
    """yield vertex lists of all cycles (each once, min vertex first)"""
    for r in range(n):
        stack = [(r, [r], {r})]
        while stack:
            v, path, used = stack.pop()
            for w in adj[v]:
                if w == r and len(path) >= 3:
                    if path[1] < path[-1]:
                        yield path[:]
                elif w > r and w not in used:
                    stack.append((w, path + [w], used | {w}))


def longest_cycles(n, adj):
    best, out = 0, []
    for cyc in all_cycles(n, adj):
        if len(cyc) > best:
            best, out = len(cyc), [cyc]
        elif len(cyc) == best:
            out.append(cyc)
    return best, out


def ears(n, adj, cyc):
    """all ears of C: (u_idx, v_idx, length, frozenset(interior))"""
    pos = {v: i for i, v in enumerate(cyc)}
    onC = set(cyc)
    out = []
    # chords
    for i, u in enumerate(cyc):
        for w in adj[u]:
            if w in onC and pos[w] != (i + 1) % len(cyc) and pos[w] != (i - 1) % len(cyc):
                if pos[w] > i:
                    out.append((i, pos[w], 1, frozenset()))
    # ears through vertices off C
    off = [v for v in range(n) if v not in onC]
    for i, u in enumerate(cyc):
        # DFS from u into the off-C region
        stack = [(w, [w], {w}) for w in adj[u] if w not in onC]
        while stack:
            v, path, used = stack.pop()
            for w in adj[v]:
                if w in onC:
                    j = pos[w]
                    if j > i:
                        out.append((i, j, len(path) + 1, frozenset(path)))
                elif w not in used:
                    stack.append((w, path + [w], used | {w}))
    # dedupe
    return sorted(set(out))


def arcs(c, i, j):
    """the two arcs between positions i<j: as (length, set of positions strictly inside)"""
    inner = set(range(i + 1, j))
    outer = set(range(j + 1, c)) | set(range(0, i))
    return (j - i, inner), (c - (j - i), outer)


def analyse(g6, verbose=True):
    n, e = from_graph6(g6)
    adj = adjacency(n, e)
    c, cycs = longest_cycles(n, adj)
    report = []
    for cyc in cycs:
        E = ears(n, adj, cyc)
        # exchanges: (delta, replaced-arc positions, interior vertices, endpoints)
        ex = []
        for (i, j, L, interior) in E:
            for (alen, inner) in arcs(c, i, j):
                d = alen - L
                if d >= 0:
                    ex.append((d, frozenset(inner), interior, (i, j)))
        small = [x for x in ex if x[0] <= 2 and x[0] > 0]
        # best independent family among small-decrement exchanges (exact, small sets)
        best = 0
        m = len(small)
        for r in range(min(m, 8), 0, -1):
            if best >= (c + 1) // 2:
                break
            for comb in itertools.combinations(range(m), r):
                ok = True
                arcsU, intU = set(), set()
                tot = 0
                for t in comb:
                    d, inner, interior, _ = small[t]
                    if (arcsU & inner) or (intU & interior):
                        ok = False; break
                    arcsU |= inner; intU |= interior; tot += d
                if ok and tot > best:
                    best = tot
        report.append((c, len(E), len(small), sorted(set(x[0] for x in ex)), best))
    c_, nE, nS, deltas, best = max(report, key=lambda r: r[4])
    if verbose:
        print("%s  n=%d  circumference=%d" % (g6, n, c_))
        print("   ears on a longest cycle: %d ;  available decrements: %s" % (nE, deltas))
        print("   exchanges with delta in {1,2}: %d ;  best independent sum: %d  (need >= %d)"
              % (nS, best, (c_ + 1) // 2))
        print("   CONJECTURE (all delta<=2, sum>=c/2): %s" % ("HOLDS" if best >= (c_ + 1) // 2 else "FAILS"))
    return best >= (c_ + 1) // 2


if __name__ == '__main__':
    if len(sys.argv) > 1:
        for g6 in sys.argv[1:]:
            analyse(g6)
    else:
        tot = fail = 0
        for line in sys.stdin:
            g6 = line.strip()
            if not g6:
                continue
            tot += 1
            if not analyse(g6, verbose=False):
                fail += 1
                if fail <= 5:
                    analyse(g6)
        print("TOTAL=%d FAILURES=%d" % (tot, fail))
