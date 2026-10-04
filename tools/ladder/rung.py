"""Rung-lattice sweep: for each escape-orbit and each leak-pair,
enumerate all K3,3 routes between the two corridors' ends (avoiding both),
collect the banned diagonal values on (u+v) and (u-v) mod 8, and report
whether the 8x8 residue lattice of rung positions dies entirely.

Lattice dead  <=>  all 8 values of u+v banned  OR  all 8 of u-v banned.
"""
from itertools import product
from sweep3p import ROWS, COLS, PERMS, ok, canon

LEAKS = [((1, 0), (0, 1)),   # (g1, g2)
         ((0, 2), (2, 0)),   # (qa, qb)
         ((1, 2), (2, 1))]   # (q'a, q'b)

def paths(M, banned, start, goal):
    """All simple paths in K3,3 (vertices ('r',i)/('c',j)) minus banned
    edges, from start to goal. Yields residue sums."""
    out = []
    def nbrs(x):
        kind, idx = x
        if kind == 'r':
            for j in COLS:
                if (idx, j) not in banned:
                    yield ('c', j), M[idx][j]
        else:
            for i in ROWS:
                if (i, idx) not in banned:
                    yield ('r', i), M[i][idx]
    def dfs(x, seen, acc):
        if x == goal:
            out.append(acc % 8)
            return
        for y, w in nbrs(x):
            if y not in seen:
                dfs(y, seen | {y}, acc + w)
    dfs(start, {start}, 0)
    return out

def lattice(M, cA, cB):
    """Return (#banned NE diagonals, #banned SE diagonals, survivors)."""
    (iA, jA), (iB, jB) = cA, cB
    gA, gB = M[iA][jA] % 8, M[iB][jB] % 8
    banned = {cA, cB}
    endsA = [(('r', iA), 'u'), (('c', jA), 'gu')]   # enter at u-side / far
    endsB = [(('r', iB), 'v'), (('c', jB), 'gv')]
    NE, SE = set(), set()
    for (eB, tagB), (eA, tagA) in product(endsB, endsA):
        for R in set(paths(M, banned, eB, eA)):
            if tagB == 'v' and tagA == 'u':
                NE.add((-1 - R) % 8)
            elif tagB == 'v' and tagA == 'gu':
                SE.add((-1 - R - gA) % 8)      # v - u ≡ c → u−v ≡ −c line
            elif tagB == 'gv' and tagA == 'u':
                SE.add((1 + R + gB) % 8)       # u − v ≡ ...
            else:
                NE.add((1 + R + gA + gB) % 8)
    cov = 64 - (8 - len(NE)) * (8 - len(SE))
    return len(NE), len(SE), 64 - cov

def main():
    sols = []
    for a, ap, t, g1, qa, qpa in product(range(8), repeat=6):
        g2 = (4 - a - ap - g1) % 8
        qb = (4 - a - qa - t) % 8
        qpb = (4 - ap - qpa - t) % 8
        M = [[a, g2, qa], [g1, ap, qpa], [qb, qpb, t]]
        if ok(M):
            sols.append(M)
    orbits = {}
    for M in sols:
        orbits.setdefault(canon(M), M)
    # escape orbits = those with some non-firing disjoint pair = no hot
    # non-anchor quad (conservation law): check quads directly
    from sweep3p import QUADPAIRS, ANCH, quadsum
    esc = []
    for c, M in orbits.items():
        hot = any(quadsum(M, rp, cp) == 4
                  for rp, cp in QUADPAIRS if (rp, cp) not in ANCH)
        if not hot:
            esc.append(M)
    print("escape orbits:", len(esc))
    hist = {}
    alive_orbits = 0
    worst = None
    for M in esc:
        surv_total = 0
        per = []
        for cA, cB in LEAKS:
            ne, se, surv = lattice(M, cA, cB)
            per.append(surv)
            surv_total += surv
        key = tuple(sorted(per))
        hist[key] = hist.get(key, 0) + 1
        if surv_total > 0:
            alive_orbits += 1
            if worst is None or surv_total > worst[0]:
                worst = (surv_total, M, per)
    print("per-orbit (sorted per-leak survivor counts) histogram:")
    for k in sorted(hist):
        print("  %s : %d orbits" % (k, hist[k]))
    print("orbits with ANY surviving rung-residue:", alive_orbits, "of", len(esc))
    if worst:
        print("worst orbit (most survivors):", worst[1], "per-leak:", worst[2])

if __name__ == "__main__":
    main()
