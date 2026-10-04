# -*- coding: utf-8 -*-
"""518: the band theorem, and the difference weapon.

(A) THE TRACE TEST IN CLOSED FORM.  For cubic G, the non-backtracking
matrix B has spectrum = {roots of x^2 - mu x + 2 : mu in spec A}
together with n/2 copies of each of +1, -1.  For even m,

    2 m c_m = tr B^m = SUM_i p_m(mu_i) + n ,

where p_m(mu) = lam^m + lam'^m for lam lam' = 2, lam + lam' = mu.
  - mu = 3 gives lam = 2, lam' = 1, so p_m(3) = 2^m + 1;
  - |mu| > 2 sqrt 2 gives lam, lam' real of the same sign, so
    p_m(mu) > 0 for even m;
  - |mu| <= 2 sqrt 2 gives |lam| = sqrt 2 and
    p_m(mu) = 2 . 2^{m/2} cos(m theta), so |p_m(mu)| <= 2 . 2^{m/2}.
Hence 2 m c_m >= 2^m + 1 + n - 2 n . 2^{m/2}, and

    2^{m/2} > 2n   ==>   c_m > 0                          (the test)

with NO Ramanujan hypothesis.  Applied with m = 2^k the unique power
of two in [g, 2g), this is the BAND THEOREM.

(B) THE DIFFERENCE WEAPON.  In any graph, two cycles meeting in
exactly one path of length p have a symmetric difference which is a
single cycle, of length a + b - 2p.  So "no 8-cycle" forbids every
pair with a + b - 2p = 8.  At a = b = g that is p = g - 4:
    g = 5 -> p = 1 : pentagons pairwise NOT sharing exactly an edge
    g = 6 -> p = 2 : hexagons  pairwise NOT sharing exactly a 2-path
    g = 7 -> p = 3 : heptagons pairwise NOT sharing exactly a 3-path
which caps c_g by a count of the sub-paths.
"""

import sys
from itertools import combinations


# ---------------- (A) the band theorem ------------------------------

def target(g):
    k = 1
    while (1 << k) < g:
        k += 1
    return k, 1 << k


def moore(g):
    if g % 2:
        return 3 * (1 << ((g - 1) // 2)) - 2
    return (1 << (g // 2 + 1)) - 2


def band_table():
    print("band theorem: girth g, band k with 2^{k-1} < g <= 2^k,")
    print("target 2^k, test fires iff 2^{2^{k-1}} > 2n i.e. "
          "n < 2^{2^{k-1}-1}")
    print()
    print(" g   k  target   Moore n>=   band bound n>=   which is "
          "stronger")
    for g in list(range(5, 20)) + [24, 26, 30, 32, 33]:
        k, t = target(g)
        thr = 1 << ((1 << (k - 1)) - 1)
        mo = moore(g)
        if g == t:
            note = "target = girth, trivial"
        elif thr > mo:
            note = "BAND (x%.1f)" % (float(thr) / mo)
        else:
            note = "Moore"
        print("%2d  %2d  %6d   %10d   %14s   %s"
              % (g, k, t, mo,
                 ("%d" % thr) if thr < 10 ** 12 else "%.3g" % thr,
                 note))
    sys.stdout.flush()


def check_452_sweep():
    """452.4 records the (n,g) test firing at n = 58, 70, 90, 102, 126
    and at 1022, 2046, 4094, 8190, 16382, 262142, 1048574, and FAILING
    at 258, 510, 65534.  Those n are Moore-extremal, so recover g from
    n and check the closed form against the record."""
    print()
    print("cross-check against the sweep recorded in 452.4")
    girth_of = {}
    for g in range(3, 70):
        girth_of[moore(g)] = g
    fires = [58, 70, 90, 102, 126, 1022, 2046, 4094, 8190, 16382,
             262142, 1048574]
    fails = [258, 510, 65534]
    for lab, lst, want in (("FIRES", fires, True),
                           ("FAILS", fails, False)):
        for n in lst:
            g = girth_of.get(n)
            if g is None:
                print("   n=%-8d  (not Moore-extremal; skipped)" % n)
                continue
            k, t = target(g)
            thr = 1 << ((1 << (k - 1)) - 1)
            got = n < thr
            trivial = (g == t)
            ok = (got == want) or (trivial and want)
            print("   n=%-8d g=%2d target=%3d  n<2^{2^{k-1}-1}=%s  "
                  "record says %s  %s%s"
                  % (n, g, t, "yes" if got else "no ", lab,
                     "AGREES" if ok else "*** DISAGREES",
                     "  (target=girth)" if trivial else ""))
    sys.stdout.flush()


# ---------------- (B) the difference weapon -------------------------

def all_cycles(adj):
    out = []
    vs = sorted(adj)
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        stack = [(s, (s,), 1 << s)]
        while stack:
            v, path, msk = stack.pop()
            for w in adj[v]:
                if w == s:
                    if len(path) >= 3 and rank[path[1]] < rank[path[-1]]:
                        out.append(path)
                elif rank[w] > rank[s] and not (msk >> w & 1):
                    stack.append((w, path + (w,), msk | (1 << w)))
    return out


def edgeset(path):
    k = len(path)
    return set(frozenset((path[i], path[(i + 1) % k]))
               for i in range(k))


def shared_is_single_path(E1, E2, V1, V2):
    """the intersection of the two cycles is exactly one path; return
    its length in edges, or None."""
    E = E1 & E2
    V = V1 & V2
    if not V:
        return None
    # the shared vertex set must induce a path in the shared edges
    deg = {}
    for e in E:
        for x in e:
            deg[x] = deg.get(x, 0) + 1
    if any(d > 2 for d in deg.values()):
        return None
    if len(V) != len(E) + 1:
        return None          # not a single path on the shared vertices
    if any(v not in deg for v in V) and len(E):
        return None
    if not E:
        return 0 if len(V) == 1 else None
    # connectivity of the shared edges
    start = next(iter(V))
    seen, dq = set([start]), [start]
    nb = {}
    for e in E:
        a, b = tuple(e)
        nb.setdefault(a, []).append(b)
        nb.setdefault(b, []).append(a)
    while dq:
        x = dq.pop()
        for y in nb.get(x, []):
            if y not in seen:
                seen.add(y)
                dq.append(y)
    if seen != V:
        return None
    return len(E)


def weapon_check(name, adj):
    cyc = all_cycles(adj)
    lens = set(len(c) for c in cyc)
    n = len(adj)
    g = min(lens)
    data = [(len(c), edgeset(c), set(c)) for c in cyc]
    bad = 0
    tested = 0
    for (a, E1, V1), (b, E2, V2) in combinations(data, 2):
        p = shared_is_single_path(E1, E2, V1, V2)
        if p is None or p == 0:
            continue
        tested += 1
        if (a + b - 2 * p) not in lens:
            bad += 1
    # the c_g cap implied by "no 8-cycle", reported for information
    cg = sum(1 for c in cyc if len(c) == g)
    cap = {5: 3.0 * n / 10, 6: n / 2.0, 7: 6.0 * n / 7}.get(g)
    print("%-16s n=%3d g=%d  |lengths|=%2d  pairs meeting in one path "
          "tested=%6d  difference-not-a-cycle-length=%d  c_g=%3d  "
          "no-8 cap=%s  8-cycle forced=%s"
          % (name, n, g, len(lens), tested, bad, cg,
             ("%.1f" % cap) if cap else "-",
             ("YES" if cap and cg > cap else "no")))
    sys.stdout.flush()
    return bad


if __name__ == '__main__':
    band_table()
    check_452_sweep()
    print()
    from cage516 import lcf, petersen, coxeter
    tot = 0
    tot += weapon_check("Petersen", petersen())
    tot += weapon_check("dodecahedron",
                        lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2))
    tot += weapon_check("Heawood", lcf([5, -5], 7))
    tot += weapon_check("Mobius-Kantor", lcf([5, -5], 8))
    tot += weapon_check("Pappus", lcf([5, 7, -7, 7, -7, -5], 3))
    print("total failures of the difference weapon:", tot)
