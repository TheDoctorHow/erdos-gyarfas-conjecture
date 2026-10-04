"""Sweep stages 2-4: hatch-cap fires, banned span-sum tables + odd-span
analogue (single-residue SAT), escaping-pair structure, corridor parities."""
from itertools import product
from sweep3p import ROWS, COLS, PERMS, ok, canon, quadsum, iw_fires

def edges():
    return [(i, j) for i in ROWS for j in COLS]

def hatch_fires(M):
    """Adjacent pairs (share row or col). Host c=(i,j), partner c'=(i,l):
    around-route from shared row-vertex i to c's far end col-j avoiding
    c, c': i -> l'' -> k' -> j with l'' != j,l ; k' != i.
    Fires iff 4 + M[i][j] + R3 == 0 (mod 8) for some route (either host)."""
    fires = 0
    pairs = []
    for i in ROWS:
        cols2 = [(j, l) for j in COLS for l in COLS if j < l]
        for j, l in cols2:
            pairs.append((('r', i), (i, j), (i, l)))
    for j in COLS:
        for i, k in [(a, b) for a in ROWS for b in ROWS if a < b]:
            pairs.append((('c', j), (i, j), (k, j)))
    for share, c1, c2 in pairs:
        hit = False
        for host, other in ((c1, c2), (c2, c1)):
            hi, hj = host
            if share[0] == 'r':
                lps = [x for x in COLS if x != hj and x != other[1]]
                for lp in lps:
                    for kp in [x for x in ROWS if x != hi]:
                        R3 = (M[hi][lp] + M[kp][lp] + M[kp][hj]) % 8
                        if (4 + M[hi][hj] + R3) % 8 == 0:
                            hit = True
            else:
                kps = [x for x in ROWS if x != hi and x != other[0]]
                for kp in kps:
                    for lp in [x for x in COLS if x != hj]:
                        R3 = (M[kp][hj] + M[kp][lp] + M[hi][lp]) % 8
                        if (4 + M[hi][hj] + R3) % 8 == 0:
                            hit = True
        if hit:
            fires += 1
    return fires, len(pairs)

def banned_dsum(M):
    """Per unordered corridor-pair: set of banned (span1+span2) mod 8."""
    B = {}
    E = edges()
    for x in range(len(E)):
        for y in range(x + 1, len(E)):
            (i, j), (k, l) = E[x], E[y]
            bans = set()
            if i != k and j != l:                     # disjoint
                bans.add((2 + quadsum(M, (i, k), (j, l))) % 8)
                k2 = [r for r in ROWS if r not in (i, k)][0]
                l2 = [c for c in COLS if c not in (j, l)][0]
                Rvv = M[k2][j] + M[k2][l]
                Ruu = M[i][l2] + M[k][l2]
                bans.add((2 + M[i][j] + M[k][l] + Rvv + Ruu) % 8)
            elif i == k:                              # share row i
                for kp in [r for r in ROWS if r != i]:
                    R2 = M[kp][l] + M[kp][j]
                    bans.add((2 + M[i][j] + M[i][l] + R2) % 8)
            else:                                     # share col j
                for lp in [c for c in COLS if c != j]:
                    R2 = M[i][lp] + M[k][lp]
                    bans.add((2 + M[i][j] + M[k][j] + R2) % 8)
            B[(E[x], E[y])] = bans
    return B

def odd_span_sat(M):
    """Single-residue-per-corridor assignment r_c with:
    2*r_c mod 8 != 6 (within-corridor disjoint table), and
    r_c + r_c' not banned for every pair. Returns #assignments (0 = dead)."""
    B = banned_dsum(M)
    E = edges()
    cand = {e: [r for r in range(8) if (2 * r) % 8 != 6] for e in E}
    count = 0
    def bt(idx, assign):
        nonlocal count
        if count > 0:
            return              # existence only
        if idx == len(E):
            count += 1
            return
        e = E[idx]
        for r in cand[e]:
            good = True
            for f, rf in assign.items():
                key = (e, f) if (e, f) in B else (f, e)
                if (r + rf) % 8 in B[key]:
                    good = False
                    break
            if good:
                assign[e] = r
                bt(idx + 1, assign)
                del assign[e]
    bt(0, {})
    return count

def escape_pairs(M):
    """Which disjoint pairs do NOT fire the inevitable web."""
    esc = []
    for i, j in product(ROWS, COLS):
        for k, l in product(ROWS, COLS):
            if k == i or l == j or (i, j) > (k, l):
                continue
            hit = False
            for host, d in ((((i, j)), (k, l)), ((k, l), (i, j))):
                hi, hj = host
                for lp in [c for c in COLS if c != hj]:
                    for kp in [r for r in ROWS if r != hi]:
                        if d in {(hi, lp), (kp, lp), (kp, hj)}:
                            continue
                        R3 = (M[hi][lp] + M[kp][lp] + M[kp][hj]) % 8
                        if (4 + M[hi][hj] + R3) % 8 == 0:
                            hit = True
            if not hit:
                esc.append(((i, j), (k, l)))
    return esc

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
    print("orbits:", len(orbits))

    hh = {}
    sat_alive = 0
    esc_kinds = {}
    evens_hist = {}
    both_full = 0
    for c, M in orbits.items():
        hf, htot = hatch_fires(M)
        hh[hf] = hh.get(hf, 0) + 1
        s = odd_span_sat(M)
        if s:
            sat_alive += 1
        ev = sum(1 for i in ROWS for j in COLS if M[i][j] % 2 == 0)
        evens_hist[ev] = evens_hist.get(ev, 0) + 1
        iwf, iwtot = iw_fires(M)
        if iwf == iwtot and hf == htot:
            both_full += 1
        if iwf < iwtot:
            esc = escape_pairs(M)
            anchors = [((0,0),(1,1)), ((0,0),(2,2)), ((1,1),(2,2))]
            kind = tuple(sorted(
                'diag' if tuple(sorted(p)) in [tuple(sorted(a)) for a in anchors]
                else 'off' for p in esc))
            esc_kinds[kind] = esc_kinds.get(kind, 0) + 1
    print("hatch fires over %d adjacent pairs, histogram:" % htot)
    for f in sorted(hh):
        print("  fires=%d : %d orbits" % (f, hh[f]))
    print("orbits with BOTH tests all-fire:", both_full)
    print("odd-span single-residue SAT alive:", sat_alive, "of", len(orbits))
    print("even-residue corridor count histogram:", dict(sorted(evens_hist.items())))
    print("escaping-pair kinds (diag = anchor-arc pair):",
          dict(sorted(esc_kinds.items())))

if __name__ == "__main__":
    main()
