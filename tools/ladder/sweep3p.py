"""3-petal CSP sweep: enumerate residue patterns, quotient by symmetry,
run Tier-2 kill recomputations (inevitable-web fire test, banned span-sum
tables) per pattern.

Matrix M[i][j], rows A={x1,x3,s2}, cols B={x2,x4,s1}:
  M = [[a, g2, qa], [g1, ap, qpa], [qb, qpb, t]]
Anchors (== 4 mod 8): Q(01|01)=Z*, Q(02|02)=Z, Q(12|12)=Z'.
All 9 quads and 6 hexagons != 0 mod 8.
"""
from itertools import permutations, product

ROWS = COLS = (0, 1, 2)
QUADPAIRS = [((0, 1), (0, 1)), ((0, 2), (0, 2)), ((1, 2), (1, 2)),
             ((0, 1), (0, 2)), ((0, 1), (1, 2)), ((0, 2), (0, 1)),
             ((0, 2), (1, 2)), ((1, 2), (0, 1)), ((1, 2), (0, 2))]
ANCH = {((0, 1), (0, 1)), ((0, 2), (0, 2)), ((1, 2), (1, 2))}
PERMS = list(permutations((0, 1, 2)))

def quadsum(M, rp, cp):
    (i, k), (j, l) = rp, cp
    return (M[i][j] + M[i][l] + M[k][j] + M[k][l]) % 8

def ok(M):
    for rp, cp in QUADPAIRS:
        s = quadsum(M, rp, cp)
        if ((rp, cp) in ANCH):
            if s != 4:
                return False
        elif s == 0:
            return False
    T = sum(sum(r) for r in M) % 8
    for p in PERMS:
        d = (M[0][p[0]] + M[1][p[1]] + M[2][p[2]]) % 8
        if (T - d) % 8 == 0:
            return False
    return True

def canon(M):
    best = None
    for p in PERMS:
        for tr in (0, 1):
            N = [[(M[p[i]][p[j]] if not tr else M[p[j]][p[i]]) for j in COLS]
                 for i in ROWS]
            t = tuple(tuple(r) for r in N)
            if best is None or t < best:
                best = t
    return best

def iw_fires(M):
    """Inevitable-web fire test per disjoint corridor pair (host, partner).
    Host c=(i,j): around-routes i->l'->k'->j (3 corridors), avoiding c and d.
    Fires iff 4 + M[i][j] + R3 == 0 (mod 8) for some valid route."""
    fires = 0
    tot = 0
    for i, j in product(ROWS, COLS):
        for k, l in product(ROWS, COLS):
            if k == i or l == j:
                continue
            if (i, j) > (k, l):
                continue
            tot += 1
            hit = False
            for host, d in (((i, j), (k, l)), ((k, l), (i, j))):
                hi, hj = host
                for lp in COLS:
                    if lp == hj:
                        continue
                    for kp in ROWS:
                        if kp == hi:
                            continue
                        route = {(hi, lp), (kp, lp), (kp, hj)}
                        if d in route:
                            continue
                        R3 = (M[hi][lp] + M[kp][lp] + M[kp][hj]) % 8
                        if (4 + M[hi][hj] + R3) % 8 == 0:
                            hit = True
                if hit:
                    break
            if hit:
                fires += 1
    return fires, tot

def main():
    sols = []
    for a, ap, t, g1, qa, qpa in product(range(8), repeat=6):
        g2 = (4 - a - ap - g1) % 8
        qb = (4 - a - qa - t) % 8
        qpb = (4 - ap - qpa - t) % 8
        M = [[a, g2, qa], [g1, ap, qpa], [qb, qpb, t]]
        if ok(M):
            sols.append(M)
    print("raw solutions:", len(sols))
    orbits = {}
    for M in sols:
        orbits.setdefault(canon(M), M)
    print("orbits (S3 x transpose):", len(orbits))
    hist = {}
    full = []
    for c, M in orbits.items():
        f, tot = iw_fires(M)
        hist[f] = hist.get(f, 0) + 1
        if f == tot:
            full.append(M)
    print("inevitable-web fires over %d disjoint pairs, histogram:" % tot)
    for f in sorted(hist):
        print("  fires=%d : %d orbits" % (f, hist[f]))
    print("orbits with ALL pairs firing (full class-caps):", len(full))
    least = min(hist)
    weak = [M for c, M in orbits.items() if iw_fires(M)[0] == least]
    print("example weakest pattern (fires=%d):" % least, weak[0] if weak else None)

if __name__ == "__main__":
    main()
