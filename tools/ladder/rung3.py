"""Stratum 3: strata-tower test. Odd-rung webs (1,3,5,7 rungs, one route)
ban full parity-classes of each stratum-1 diagonal line; even-rung webs
(2,4,6,8 rungs, route-pairs) kill the (0,0) class iff the pair-menu meets
parity gA+gB. A leak survives all strata iff:
  (a) every pair-route sum has parity != (gA+gB) mod 2, AND
  (b) the odd-union leaves a nonempty class grid:
      all NE single-route lines share one parity and all SE lines share one
      parity (survivor classes = complementary parities, consistent grid).
"""
from itertools import product
from sweep3p import ROWS, COLS, ok, canon, QUADPAIRS, ANCH, quadsum
from rung2 import all_paths, pair_menu, LEAKS

def single_route_lines(M, cA, cB):
    """Stratum-1 diagonal values: NE (u+v) and SE (u-v) line constants."""
    (iA, jA), (iB, jB) = cA, cB
    gA, gB = M[iA][jA] % 8, M[iB][jB] % 8
    banned = {cA, cB}
    Ar, Ac = ('r', iA), ('c', jA)
    Br, Bc = ('r', iB), ('c', jB)
    NE, SE = set(), set()
    for R, _ in all_paths(M, banned, Br, Ar):
        NE.add((-1 - R) % 8)
    for R, _ in all_paths(M, banned, Br, Ac):
        SE.add((-1 - R - gA) % 8)
    for R, _ in all_paths(M, banned, Bc, Ar):
        SE.add((1 + R + gB) % 8)
    for R, _ in all_paths(M, banned, Bc, Ac):
        NE.add((1 + R + gA + gB) % 8)
    return NE, SE

def leak_survives(M, cA, cB):
    gA = M[cA[0]][cA[1]] % 8
    gB = M[cB[0]][cB[1]] % 8
    # even-strata: fires iff pair-menu contains parity gA+gB
    pm = pair_menu(M, cA, cB)
    if any(x % 2 == (gA + gB) % 2 for x in pm):
        return False
    # odd-strata union: parity classes of single-route lines
    NE, SE = single_route_lines(M, cA, cB)
    ne_par = {x % 2 for x in NE}
    se_par = {x % 2 for x in SE}
    if len(ne_par) == 2 or len(se_par) == 2:
        return False              # odd-union covers everything
    # survivor classes: u+v parity != ne_par, u-v parity != se_par;
    # grid consistency: (u+v)+(u-v) even
    pne = 1 - next(iter(ne_par)) if ne_par else None
    pse = 1 - next(iter(se_par)) if se_par else None
    if pne is not None and pse is not None and (pne + pse) % 2 != 0:
        return False              # inconsistent grid: no classes
    return True

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
    esc = [M for c, M in orbits.items()
           if not any(quadsum(M, rp, cp) == 4
                      for rp, cp in QUADPAIRS if (rp, cp) not in ANCH)]
    hist = {}
    alive = []
    for M in esc:
        n = sum(1 for cA, cB in LEAKS if leak_survives(M, cA, cB))
        hist[n] = hist.get(n, 0) + 1
        if n > 0:
            alive.append((M, n))
    print("surviving leaks per orbit, histogram:", dict(sorted(hist.items())))
    print("orbits with ANY leak surviving all strata:", len(alive),
          "of", len(esc))
    if alive:
        worst = max(alive, key=lambda x: x[1])
        print("example:", worst[0], "surviving leaks:", worst[1])

if __name__ == "__main__":
    main()
