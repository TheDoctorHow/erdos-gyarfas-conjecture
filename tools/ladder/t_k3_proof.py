"""MACHINE-CHECKED PROOF that T(2N) has no three-chord generator for
all large N (N = 0 mod 3).

T(2N) = { S_j = (3j, 3j+4) } union { D_j = (3j+2, 3j+2+N) }.

A three-chord generator instance (597.2) is: three disjoint sites with
gap profile a placement of (3,1,1) or (2,2,1), boundaries wired by one
of the four valid pairings, every chord a chord of T.  Each valid
pairing wires the three sites as a triangle - one chord per site pair.

Positions are represented SYMBOLICALLY as c + e*N (c a small integer,
e an integer reduced mod 2, since 2N = 0).  Every chord of T imposes an
exact offset between its endpoints:

    S:  higher = lower + 4,  lower = 0 mod 3
    D:  other  = one  + N,   one   = 2 mod 3   (then other = 2 too)

Fixing site 1's left boundary at 0 with unknown residue r in {0,1,2},
two chords of the triangle determine sites 2 and 3; the third chord
must close: an equation  c + e*N = 0 mod 2N.  For N > 12 this forces
e even and c = 0 exactly.  Residues are checked mod 3 (N = 0 mod 3, so
residue(c + eN) = residue(c) + r-part).  Surviving templates are then
checked for site disjointness: two sites with equal e-part must have
disjoint small-integer intervals; unequal e-parts are N apart, disjoint
for large N.

If NO template survives, T(2N) has no three-chord generator for every
N beyond the small-N threshold - a finite machine-checked case
analysis, 6 profiles x 4 pairings x 27 chord-type assignments.
"""
import sys
from itertools import product, permutations

PAIRINGS = [((0,2),(1,4),(3,5)), ((0,3),(1,4),(2,5)),
            ((0,3),(1,5),(2,4)), ((0,4),(1,3),(2,5))]
PROFILES = sorted(set(permutations((3,1,1))) | set(permutations((2,2,1))))

def site_of(b): return b // 2
def side_of(b): return b % 2          # 0 = left, 1 = right

def check_template(prof, pairing, types):
    """types: per chord, one of 'S01' (b = a+4, a class 0),
    'S10' (a = b+4, b class 0), 'D' (b = a+N, a class 2; symmetric).
    Chord i joins boundary pair pairing[i] = (a_b, b_b).
    Returns a witness dict if the template is realisable for all large
    N, else None."""
    # unknown left-boundary positions of sites: P[s] = (c, e); P[0] fixed
    # unknown residue r of P[0] mod 3; other residues derived
    for r0 in range(3):
        P = {0: (0, 0)}
        R = {0: r0}                      # residue of site s's left boundary
        # process chords as constraints; iterate to propagate (triangle:
        # two chords fix sites 2,3, third closes)
        cons = []
        ok = True
        for (ab, bb), t in zip(pairing, types):
            cons.append((ab, bb, t))
        changed = True
        closure_ok = True
        guard = 0
        pend = list(cons)
        while pend and guard < 10:
            guard += 1
            nxt = []
            for ab, bb, t in pend:
                sa, sb = site_of(ab), site_of(bb)
                ga = prof[sa] if side_of(ab) else 0
                gb = prof[sb] if side_of(bb) else 0
                # position of boundary = P[site] + (gap if right side)
                known_a, known_b = sa in P, sb in P
                if not (known_a or known_b):
                    nxt.append((ab, bb, t)); continue
                # offset o with pos_b = pos_a + o  (as (c, e))
                if t == 'S01':   o = (4, 0)
                elif t == 'S10': o = (-4, 0)
                else:            o = (0, 1)
                if known_a and not known_b:
                    ca, ea = P[sa]
                    P[sb] = (ca + ga + o[0] - gb, ea + o[1])
                elif known_b and not known_a:
                    cb, eb = P[sb]
                    P[sa] = (cb + gb - o[0] - ga, eb - o[1])
                else:
                    ca, ea = P[sa]; cb, eb = P[sb]
                    dc = (cb + gb) - (ca + ga) - o[0]
                    de = eb - ea - o[1]
                    # closure: dc + de*N = 0 mod 2N, N large
                    # => de even and dc = 0   (de odd needs dc = ±N: fails)
                    if not (dc == 0 and de % 2 == 0):
                        closure_ok = False
            pend = nxt
        if pend or not closure_ok:
            continue
        # residues: residue of P[s] left boundary = (c + 0) mod 3 shifted
        # by r0: define res(s) = (r0 + P[s].c - P[0].c) mod 3  (N = 0 mod 3)
        def res_boundary(b):
            s = site_of(b)
            g = prof[s] if side_of(b) else 0
            return (r0 + P[s][0] + g) % 3
        res_ok = True
        for (ab, bb), t in zip(pairing, types):
            if t == 'S01':
                if res_boundary(ab) != 0: res_ok = False
            elif t == 'S10':
                if res_boundary(bb) != 0: res_ok = False
            else:
                if res_boundary(ab) != 2: res_ok = False
        if not res_ok:
            continue
        # site disjointness for large N: same e-part => integer intervals
        # must be disjoint; different e-part (mod 2) => N apart, fine
        disj = True
        sites = []
        for s in range(3):
            c, e = P[s]
            sites.append((c, e % 2, prof[s]))
        for i in range(3):
            for j in range(i + 1, 3):
                ci, ei, gi = sites[i]
                cj, ej, gj = sites[j]
                if ei == ej:
                    if not (ci + gi < cj or cj + gj < ci):
                        disj = False
        if not disj:
            continue
        return {'r0': r0, 'sites': sites}
    return None

if __name__ == '__main__':
    survivors = []
    total = 0
    for prof in PROFILES:
        for pairing in PAIRINGS:
            for types in product(('S01', 'S10', 'D'), repeat=3):
                total += 1
                w = check_template(prof, pairing, types)
                if w:
                    survivors.append((prof, pairing, types, w))
    print('templates checked:', total)
    print('templates realisable for all large N:', len(survivors))
    for s in survivors:
        print('  SURVIVOR:', s)
    if not survivors:
        print()
        print('PROVED: T(2N) has no three-chord generator for all')
        print('sufficiently large N with N = 0 mod 3 (every template')
        print('fails closure, residue, or disjointness).')
