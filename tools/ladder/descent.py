"""LOCAL DESCENT LEMMA test.

The Ear Subset-Sum Conjecture failed because it fixed ONE base cycle (a longest one)
and asked for many small independent decrements.  Decrements off a longest cycle are
few and large.  The repair: don't fix the base.  Descend one step at a time, using the
ears of the CURRENT cycle -- a shorter cycle has more vertices outside it, hence more
bridges and more ears, so descent should get EASIER as we go down.

  LOCAL DESCENT LEMMA (under test).  Let G be a 3-connected cubic graph with
  circumference c, and let C be a cycle with |C| even and |C| > c/2.  Then some single
  ear exchange on C yields a cycle of length |C| - 2.

An ear P of C with endpoints splitting C into arcs of lengths a and |C|-a produces
cycles of lengths a+|P| and (|C|-a)+|P|.  So a one-step descent by 2 exists iff some
ear has |P| = a-2 for one of its two arcs a.

Reported separately: the weaker "descend by 1 or 2", and how often descent needs a
non-local move.
"""
import sys, collections
sys.path.insert(0, '.')
from ref import from_graph6


def adjacency(n, edges):
    a = collections.defaultdict(list)
    for u, v in edges:
        a[u].append(v); a[v].append(u)
    return a


def all_cycles(n, adj):
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


def three_edge_connected(n, edges):
    m = len(edges)
    for i in range(m):
        for j in range(-1, m):
            if j != -1 and j <= i:
                continue
            rem = [e for k, e in enumerate(edges) if k != i and k != j]
            a = adjacency(n, rem)
            seen = {0}; st = [0]
            while st:
                u = st.pop()
                for w in a[u]:
                    if w not in seen:
                        seen.add(w); st.append(w)
            if len(seen) != n:
                return False
    return True


def ear_decrements(n, adj, cyc):
    """all achievable new-cycle lengths from a single ear exchange on cyc"""
    c = len(cyc)
    pos = {v: i for i, v in enumerate(cyc)}
    onC = set(cyc)
    out = set()
    # chords
    for i, u in enumerate(cyc):
        for w in adj[u]:
            if w in onC:
                j = pos[w]
                d = (j - i) % c
                if d not in (0, 1, c - 1):
                    out.add(d + 1); out.add(c - d + 1)
    # ears through off-cycle vertices
    for i, u in enumerate(cyc):
        stack = [(w, 1, {w}) for w in adj[u] if w not in onC]
        while stack:
            v, L, used = stack.pop()
            for w in adj[v]:
                if w in onC:
                    j = pos[w]
                    if j == i:
                        continue
                    d = (j - i) % c
                    out.add(d + L + 1); out.add(c - d + L + 1)
                elif w not in used:
                    stack.append((w, L + 1, used | {w}))
    return out


def analyse(g6):
    n, e = from_graph6(g6)
    if not three_edge_connected(n, e):
        return None
    adj = adjacency(n, e)
    cycles = list(all_cycles(n, adj))
    c = max(len(x) for x in cycles)
    bad2 = bad12 = tested = 0
    worst = None
    for cyc in cycles:
        L = len(cyc)
        if L % 2 or L <= c / 2:
            continue
        tested += 1
        reach = ear_decrements(n, adj, cyc)
        if L - 2 not in reach:
            bad2 += 1
            if worst is None:
                worst = (g6, n, c, L, sorted(reach))
        if L - 2 not in reach and L - 1 not in reach:
            bad12 += 1
    return tested, bad2, bad12, worst


if __name__ == '__main__':
    T = B2 = B12 = 0
    shown = 0
    for line in sys.stdin:
        g6 = line.strip()
        if not g6:
            continue
        r = analyse(g6)
        if r is None:
            continue
        t, b2, b12, worst = r
        T += t; B2 += b2; B12 += b12
        if b2 and shown < 6:
            shown += 1
            print("NO ONE-STEP -2 : %s n=%d c=%d  from cycle length %d ; reachable lengths %s"
                  % (worst[0], worst[1], worst[2], worst[3], worst[4]))
    print("even cycles above c/2 tested: %d" % T)
    print("  no single-ear descent by 2      : %d  (%.2f%%)" % (B2, 100.0 * B2 / max(T, 1)))
    print("  no single-ear descent by 1 or 2 : %d  (%.2f%%)" % (B12, 100.0 * B12 / max(T, 1)))
