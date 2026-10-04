"""hexagon438.py -- the hexagon world at m = 3, t = 3: exhaustive
feasibility of a single hexagon-bridge.

Under the cascade conditions (437.1) at m = 3, every bridge's
attachment set is one full class mod 2^(t-1): SIX positions, equally
spaced at L/6 = 2^(t-1) - a regular hexagon on C. At t = 3 (L = 24,
spacing 4, arcs in {4, 8, 12}) the shallowness table for a
within-bridge ear over arc a is:
    arc 4:  T =/= 4 (mod 8);  arc 8:  T =/= 0 (mod 8);
    arc 12: T =/= 4 (mod 8)
(both arc orientations included). Local minimality: T > 12 for every
pair - satisfiable by magnitude, so the BINDING system is mod 8.

Subcubic trees with 6 leaves have two internal topologies:
  CATERPILLAR: internal path s1-s2-s3-s4, leaves (2,1,1,2);
  SNOWFLAKE:   internal star c-(s1,s2,s3), two leaves per s_i.
Caterpillar leaf distances are SEPARABLE: T_ab = B_a + A_b for a
left of b, with A_i + B_i = 2 l_i (even). Snowflake: D_a = l_a + e_i;
T same-group = D_a + D_b - 2 e_i, cross-group = D_a + D_b.

Search: all hexagon assignments x all mod-8 value vectors, exhaustive
with pruning. Output: the number of feasible configurations.
"""
import itertools

F = {4: {4}, 8: {0}, 12: {4}}      # forbidden T mod 8 per arc
HEX = [0, 4, 8, 12, 16, 20]

def arc(p, q):
    d = abs(p - q) % 24
    return min(d, 24 - d)

def forb(p, q):
    return F[arc(p, q)]

# ---------------- caterpillar ----------------
# spine slots 0..5 (leaf order along the spine; slots 0,1 on s1... the
# separable form only needs the ORDER). assignment: slot -> hexagon pos.
def caterpillar_feasible():
    count = 0
    for perm in itertools.permutations(HEX[1:]):
        assign = (HEX[0],) + perm          # fix slot 0 -> position 0
        # DFS over slots in reverse: choose A_b for b = 5..0, and check
        # B_a feasibility once all A_b (b > a) fixed.
        # A_a, B_a in Z8 with A_a + B_a even <=> A_a == B_a (mod 2).
        sols = [dict()]
        ok_assign = False
        def dfs(a, Avals):
            # Avals: dict b -> A_b for b > a
            if a < 0:
                return True
            for Aa in range(8):
                # B_a must avoid f - A_b for all b > a, f in forb(a,b),
                # and B_a == A_a (mod 2)
                bad = set()
                for b, Ab in Avals.items():
                    for f in forb(assign[a], assign[b]):
                        bad.add((f - Ab) % 8)
                cands = [B for B in range(Aa % 2, 8, 2) if B not in bad]
                if not cands:
                    continue
                Avals2 = dict(Avals); Avals2[a] = Aa
                if dfs(a - 1, Avals2):
                    return True
            return False
        if dfs(5, {}):
            count += 1
    return count

# ---------------- snowflake ----------------
def snowflake_feasible():
    count = 0
    groups = [(0, 1), (2, 3), (4, 5)]
    for perm in itertools.permutations(HEX[1:]):
        assign = (HEX[0],) + perm
        found = False
        for Ds in itertools.product(range(8), repeat=6):
            for Es in itertools.product(range(8), repeat=3):
                ok = True
                for gi, (x, y) in enumerate(groups):
                    T = (Ds[x] + Ds[y] - 2 * Es[gi]) % 8
                    if T in forb(assign[x], assign[y]): ok = False; break
                    # parity: l_x = D_x - e_i must be a nonneg length -
                    # magnitudes free; no parity constraint here since
                    # l can be any integer >= 1
                if not ok: continue
                for (gi, (x, _)), (gj, (y, _)) in [(a, b) for a in
                        enumerate(groups) for b in enumerate(groups)
                        if a[0] < b[0]]:
                    pass
                for gi in range(3):
                    for gj in range(gi + 1, 3):
                        for x in groups[gi]:
                            for y in groups[gj]:
                                T = (Ds[x] + Ds[y]) % 8
                                if T in forb(assign[x], assign[y]):
                                    ok = False; break
                            if not ok: break
                        if not ok: break
                    if not ok: break
                if ok:
                    found = True; break
            if found: break
        count += found
    return count

c = caterpillar_feasible()
print('caterpillar hexagon-bridges: feasible assignments: %d / 120' % c)
s = snowflake_feasible()
print('snowflake   hexagon-bridges: feasible assignments: %d / 120' % s)
if c == 0 and s == 0:
    print()
    print('THE HEXAGON WORLD AT t = 3, m = 3 IS EMPTY: no subcubic')
    print('6-leaf bridge satisfies the shallowness table on a hexagon.')
else:
    print()
    print('The hexagon world is inhabited: the mod-8 system is')
    print('satisfiable; recorded honestly with a witness to be')
    print('extracted for the next stride.')
