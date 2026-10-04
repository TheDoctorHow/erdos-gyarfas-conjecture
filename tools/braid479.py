"""479 audit: Mobius-ladder normal form of a braid, and the 479.4 witness.

Three independent checks, none of them reading a printed claim back:

 (1) BRAID RIGIDITY.  Random pairwise-crossing chord families on a
     cycle are built by brute force (all pairings of 2k points that
     are pairwise crossing) and the antipodal pattern i <-> i+k is
     verified to be the ONLY one.
 (2) THE CYCLE FORMS of 479.3 (r = 0,1,2 and the separable normal
     form g_i + f_j) are checked against an exhaustive enumeration of
     the cycles of a weighted Mobius ladder.
 (3) THE 479.4 WITNESS: k = 12, all rim arcs 1, all rungs 65.  Every
     cycle of ML_12 is enumerated and its weight tested for being a
     power of two.

Cycles are enumerated from the cycle space: every element of GF(2)^E
in the span is tested for being a single cycle (connected, all
degrees 2).  That is independent of the length formulas of 479.3.
"""

from itertools import combinations


def mobius_ladder(k):
    """Vertices 0..2k-1; rim edges (j, j+1 mod 2k); rungs (i, i+k)."""
    rim = [(j, (j + 1) % (2 * k)) for j in range(2 * k)]
    rungs = [(i, i + k) for i in range(k)]
    return rim, rungs


def all_cycles(nv, edges):
    """Every cycle, as an edge bitmask.  DFS enumeration: a cycle is
    generated once, from its smallest vertex, in the direction whose
    second vertex is smaller than its last vertex."""
    eid = {}
    adj = {v: [] for v in range(nv)}
    for i, (u, v) in enumerate(edges):
        eid[(u, v)] = i
        eid[(v, u)] = i
        adj[u].append(v)
        adj[v].append(u)
    out = []

    def walk(start, path, used, mask):
        last = path[-1]
        for nxt in adj[last]:
            if nxt == start and len(path) >= 3:
                if path[1] < path[-1]:
                    out.append(mask | (1 << eid[(last, nxt)]))
            elif nxt > start and nxt not in used:
                walk(start, path + [nxt], used | {nxt},
                     mask | (1 << eid[(last, nxt)]))

    for s in range(nv):
        walk(s, [s], {s}, 0)
    return out


# ---------------------------------------------------------------- (1)

def check_rigidity(kmax=6):
    """Only the antipodal pairing is pairwise crossing."""
    bad = 0
    for k in range(2, kmax + 1):
        pts = list(range(2 * k))
        pairings = []

        def build(rem, acc):
            if not rem:
                pairings.append(tuple(acc))
                return
            a = rem[0]
            for i in range(1, len(rem)):
                b = rem[i]
                build([x for x in rem[1:] if x != b], acc + [(a, b)])

        build(pts, [])
        crossing = []
        for pr in pairings:
            ok = True
            for (a, b), (c, d) in combinations(pr, 2):
                lo, hi = min(a, b), max(a, b)
                inside_c = lo < c < hi
                inside_d = lo < d < hi
                if inside_c == inside_d:
                    ok = False
                    break
            if ok:
                crossing.append(pr)
        want = tuple(sorted((i, i + k) for i in range(k)))
        got = [tuple(sorted(tuple(sorted(e)) for e in pr))
               for pr in crossing]
        if got != [want]:
            bad += 1
            print("  k=%d RIGIDITY FAILS: %r" % (k, got))
        else:
            print("  k=%2d  pairwise-crossing pairings: %d (antipodal)"
                  % (k, len(got)))
    return bad


# ---------------------------------------------------------------- (2)

def check_forms(k, arcs, rungs):
    """479.3's r = 0,1,2 forms and separable normal form vs brute force."""
    rim, rg = mobius_ladder(k)
    edges = rim + rg
    w = arcs + rungs
    cyc = all_cycles(2 * k, edges)
    nrim = len(rim)

    def weight(mask):
        return sum(w[i] for i in range(len(edges)) if mask >> i & 1)

    def rungs_used(mask):
        return [i for i in range(k) if mask >> (nrim + i) & 1]

    L = sum(arcs)
    A = [sum(arcs[(i + s) % (2 * k)] for s in range(k)) for i in range(2 * k)]
    e = [arcs[s] + arcs[s + k] for s in range(k)]
    b = [sum(e[:i]) for i in range(k)]
    f = [rungs[i] + b[i] for i in range(k)]
    g = [rungs[i] - b[i] for i in range(k)]

    pred0 = {L}
    pred1 = set()
    for i in range(k):
        pred1.add(rungs[i] + A[i])
        pred1.add(rungs[i] + L - A[i])
    pred2 = set()
    for i in range(k):
        for j in range(i + 1, k):
            B = sum(e[i:j])
            pred2.add(rungs[i] + rungs[j] + B)
            pred2.add(rungs[i] + rungs[j] + L - B)
    sep2 = set()
    for i in range(k):
        for j in range(i + 1, k):
            sep2.add(g[i] + f[j])
            sep2.add(f[i] + g[j] + L)

    act = {0: set(), 1: set(), 2: set()}
    grading_ok = True
    for mask in cyc:
        r = len(rungs_used(mask))
        val = weight(mask)
        rho = val - sum(rungs[i] for i in rungs_used(mask))
        if not (0 <= rho <= L):
            grading_ok = False
        if r <= 2:
            act[r].add(val)

    ok = (act[0] == pred0 and act[1] == pred1 and act[2] == pred2
          and sep2 == pred2 and grading_ok)
    print("  k=%2d L=%3d cycles=%5d  r0 %s r1 %s r2 %s sep %s grade %s"
          % (k, L, len(cyc), act[0] == pred0, act[1] == pred1,
             act[2] == pred2, sep2 == pred2, grading_ok))
    return ok


# ---------------------------------------------------------------- (3)

def is_pow2(x):
    return x > 0 and x & (x - 1) == 0


def check_witness(k=12, D=65):
    rim, rg = mobius_ladder(k)
    edges = rim + rg
    w = [1] * len(rim) + [D] * k
    cyc = all_cycles(2 * k, edges)
    lens = sorted({sum(w[i] for i in range(len(edges)) if mask >> i & 1)
                   for mask in cyc})
    hits = [x for x in lens if is_pow2(x)]
    print("  witness k=%d D=%d : %d cycles, %d distinct lengths in [%d,%d]"
          % (k, D, len(cyc), len(lens), lens[0], lens[-1]))
    print("  powers of two in the spectrum: %r" % (hits,))
    return hits == []


if __name__ == "__main__":
    print("(1) BRAID RIGIDITY")
    bad = check_rigidity(6)
    print("(2) CYCLE FORMS 479.3")
    ok = True
    ok &= check_forms(4, [1, 2, 1, 3, 1, 1, 2, 1], [9, 11, 10, 12])
    ok &= check_forms(5, [1, 1, 2, 1, 1, 1, 2, 1, 1, 3],
                      [13, 12, 14, 11, 15])
    ok &= check_forms(6, [1] * 12, [17, 18, 19, 20, 21, 22])
    print("(3) THE 479.4 WITNESS")
    w = check_witness()
    print("ALL: rigidity %s, forms %s, witness %s"
          % (bad == 0, ok, w))
