"""Stratum 2: same-residue kill test. Leak dies (caps at <=64 rungs) iff
some vertex-disjoint skeleton route-pair satisfies
   R1 + R2 == -2 - gA - gB  (mod 8).
Route-pair patterns: {(Br->Ac),(Bc->Ar)} and {(Br->Ar),(Bc->Ac)},
paths in K3,3 minus the two leak edges, vertex-disjoint."""
from itertools import product
from sweep3p import ROWS, COLS, ok, canon, QUADPAIRS, ANCH, quadsum

LEAKS = [((1, 0), (0, 1)), ((0, 2), (2, 0)), ((1, 2), (2, 1))]

def all_paths(M, banned, start, goal, blocked=frozenset()):
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
    def dfs(x, seen, acc, verts):
        if x == goal:
            out.append((acc % 8, frozenset(verts)))
            return
        for y, w in nbrs(x):
            if y not in seen and y not in blocked:
                dfs(y, seen | {y}, acc + w, verts | {y})
    dfs(start, {start}, 0, {start})
    return out

def pair_menu(M, cA, cB):
    """Menu of (R1+R2) mod 8 over vertex-disjoint route pairs, both patterns."""
    (iA, jA), (iB, jB) = cA, cB
    banned = {cA, cB}
    Ar, Ac = ('r', iA), ('c', jA)
    Br, Bc = ('r', iB), ('c', jB)
    menu = set()
    for (s1, t1), (s2, t2) in [((Br, Ac), (Bc, Ar)), ((Br, Ar), (Bc, Ac))]:
        for R1, V1 in all_paths(M, banned, s1, t1):
            for R2, _ in all_paths(M, banned, s2, t2, blocked=V1):
                menu.add((R1 + R2) % 8)
    return menu

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
    print("escape orbits:", len(esc))
    hist = {}
    survivors = []
    for M in esc:
        killed = 0
        for cA, cB in LEAKS:
            gA = M[cA[0]][cA[1]] % 8
            gB = M[cB[0]][cB[1]] % 8
            target = (-2 - gA - gB) % 8
            if target in pair_menu(M, cA, cB):
                killed += 1
        hist[killed] = hist.get(killed, 0) + 1
        if killed < 3:
            survivors.append((M, killed))
    print("leaks killed per orbit, histogram:", dict(sorted(hist.items())))
    print("orbits with ALL 3 leaks capped:", hist.get(3, 0))
    print("orbits with surviving leaks:", len(survivors))
    if survivors:
        worst = min(survivors, key=lambda x: x[1])
        print("example survivor:", worst[0], "leaks killed:", worst[1])

if __name__ == "__main__":
    main()
