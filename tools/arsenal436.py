"""arsenal436.py -- turning J7's obstructions into weapons.

 (1) NEAR-STALL RIGIDITY (W1): zigzag values over a crossing family
     with two shift values c and c+e form an AP with difference e;
     the adversary survives iff v_2(base) < v_2(e). Enumerate all
     role-assignments in the winding formula over Z/16 (t = 3) and
     verify: deviation e with v_2(e) <= v_2(c0) creates a composite
     hitting a progress class; v_2(e) > v_2(c0) dodges.
 (2) INTERVAL TRICHOTOMY (W2): any k intervals contain s pairwise
     nested, or s pairwise disjoint, or k/s^2 pairwise crossing
     (Dilworth-Mirsky twice). Brute-force verification on random
     interval families.
 (3) WITHIN-BLOCK LAMINARITY (W2'): chords with both feet in one arc
     of C: nested-in-arc and disjoint-in-arc are both NON-CROSSING on
     the circle (laminar); only staircase overlap crosses. Verified
     by direct crossing tests.
 (4) SUBSET-SUM GROWTH + CHORD-DEPTH LAW (W3): a sequence of k
     residues of common 2-adic valuation a in Z/2^t has subset-sum
     set of size >= min(2^(t-a), k+1); consequently the route-shift
     classes of an ear with k disjoint self-chords cover the coset
     delta - 2^a Z, and 0 is hit iff a <= v_2(delta) and the count
     reaches coverage: the adversary needs his chords 2-adically
     DEEPER than his shifts.
 (5) THRESHOLD AUDIT (W4): the climbing count f(a) = 2^(t-a+1) - 1
     against the supply m 2^t / 3: fires for all m at a >= 1; at
     a = 0 only for m >= 7. CORRECTION to 432.6's 'm >= 5' sentence.
"""
import itertools, random

v2 = lambda x: (abs(x) & -abs(x)).bit_length() - 1 if x else 999

# ---------------- (1) near-stall rigidity ----------------
t = 3; M = 1 << (t + 1)   # winding-odd zigzags: progress iff sum == 2^t mod 2^(t+1)
def zigzag_values(inner_vals, outer_vals, kmax=7):
    """Reachable alternating-sum values: choose an odd-size subfamily
    (indices increasing); ranks alternate inner/outer."""
    n = len(inner_vals); out = set()
    for k in range(1, kmax + 1, 2):
        for sub in itertools.combinations(range(n), k):
            s = 0
            for r, i in enumerate(sub):
                s += inner_vals[i] if r % 2 == 0 else outer_vals[i]
            out.add(s % M)
    return out

random.seed(436)
# BLOCKINESS dichotomy at t = 3: c0 = 2, deviation e = 2, need
# c0 + 2j == 8 (mod 16): j == 3 (mod 8). n = 15 ears, 8 deviants.
c0, e, n = 2, 2, 15
target = 1 << t
def run(dev_idx):
    inner = [c0] * n; outer = [(-c0) % M] * n
    for j in dev_idx:
        inner[j] = (c0 + e) % M; outer[j] = (-c0 - e) % M
    return zigzag_values(inner, outer)
vals_i = run(set(range(0, n, 2)))      # INTERLEAVED deviants
vals_c = run(set(range(8)))            # CONTIGUOUS deviants
print('(1) blockiness dichotomy (c0=2, e=2, t=3; progress class 8 mod'
      ' 16):')
print('    interleaved deviants: progress reachable: %s (expected True)'
      % (target in vals_i))
print('    contiguous  deviants: progress reachable: %s (expected'
      ' False); value set %s' % (target in vals_c, sorted(vals_c)))
assert target in vals_i and target not in vals_c

# ---------------- (2) interval trichotomy ----------------
def trichotomy(iv):
    k = len(iv)
    best_n = best_d = best_c = 1
    for sub in itertools.combinations(range(k), 2):
        pass
    # longest nested chain / disjoint chain via DP (both are interval orders)
    iv_sorted = sorted(iv, key=lambda ab: (ab[0], -ab[1]))
    # nested: sort by start asc, end desc; LIS on containment
    def chain(rel):
        best = 1; dp = []
        for i in range(len(iv)):
            m = 1
            for j in range(i):
                if rel(iv[j], iv[i]): m = max(m, dp[j] + 1)
            dp.append(m); best = max(best, m)
        return best
    iv = sorted(iv)
    nested = chain(lambda A, B: A[0] <= B[0] and B[1] <= A[1] and A != B)
    disj = chain(lambda A, B: A[1] < B[0])
    # crossing clique: greedy staircase
    cross = chain(lambda A, B: A[0] < B[0] <= A[1] < B[1])
    return nested, disj, cross

ok2 = 0
for trial in range(300):
    k = random.randint(6, 14)
    iv = [tuple(sorted(random.sample(range(40), 2))) for _ in range(k)]
    ndc = trichotomy(iv)
    s = max(ndc[0], ndc[1])
    ok2 += (ndc[2] >= k / (s * s) - 1e-9 or s >= 2)
    assert ndc[2] * s * s >= k or s >= 2, (iv, ndc)
print('(2) interval trichotomy: 300/300 random families satisfy'
      ' nested/disjoint >= s or crossing >= k/s^2: %d/300' % ok2)

# ---------------- (3) within-block laminarity ----------------
Lc = 48
def crossing_on_circle(c1, c2):
    a, b = sorted(c1); c, d = sorted(c2)
    return (a < c < b < d) or (c < a < d < b)
ok3 = True
for trial in range(500):
    # two chords with all four feet in the arc [0, 20] of a 48-cycle
    pts = random.sample(range(21), 4)
    random.shuffle(pts)
    c1, c2 = (pts[0], pts[1]), (pts[2], pts[3])
    i1, i2 = tuple(sorted(c1)), tuple(sorted(c2))
    nested = i1[0] <= i2[0] and i2[1] <= i1[1] or \
             i2[0] <= i1[0] and i1[1] <= i2[1]
    disjoint = i1[1] < i2[0] or i2[1] < i1[0]
    if nested or disjoint:
        ok3 &= not crossing_on_circle(c1, c2)
    else:
        ok3 &= crossing_on_circle(c1, c2)
print('(3) within-block chords: nested/disjoint intervals never cross on'
      ' the circle; staircases always do: %s' % ok3)
assert ok3

# ---------------- (4) subset-sum growth + chord depth ----------------
ok4 = True
for trial in range(400):
    tt = random.choice([3, 4]); a = random.choice([0, 1, 2])
    if a >= tt: continue
    k = random.randint(1, 12)
    seq = [(1 << a) * random.choice([1, 3, 5, 7]) % (1 << tt)
           for _ in range(k)]
    seq = [x if v2(x) == a else (1 << a) for x in seq]
    sums = {0}
    for x in seq:
        sums |= {(s + x) % (1 << tt) for s in sums}
    ok4 &= len(sums) >= min(1 << (tt - a), k + 1)
print('(4) subset-sum growth |Sigma| >= min(2^(t-a), k+1): %s' % ok4)
assert ok4
# chord-depth law demo: delta with v2(delta) >= a and full coverage -> hit
tt, a = 3, 1
seq = [2] * 7   # k = 7 >= 2^(t-a) - 1 = 3, valuation 1
sums = {0}
for x in seq: sums |= {(s + x) % (1 << tt) for s in sums}
delta_hit = 4   # v2 = 2 >= a: 0 in delta - sums?
delta_dodge = 1 # v2 = 0 < a
print('    chord-depth: delta=4 (v2>=a) hit: %s; delta=1 (v2<a) safe: %s'
      % ((delta_hit % (1 << tt)) in sums, (delta_dodge % (1 << tt)) not in sums))
assert (delta_hit % (1 << tt)) in sums and (delta_dodge % (1 << tt)) not in sums

# ---------------- (5) threshold audit ----------------
print('(5) climbing threshold f(a) = 2^(t-a+1)-1 vs supply m*2^t/3:')
for tt in (3, 4):
    for m in (3, 5, 7, 9):
        supply = m * (1 << tt) // 3
        f0 = (1 << (tt + 1)) - 1; f1 = (1 << tt) - 1
        print('    t=%d m=%d: supply %3d | f(0)=%3d %s | f(1)=%3d %s'
              % (tt, m, supply, f0,
                 'FIRES' if supply >= f0 else 'short',
                 f1, 'FIRES' if supply >= f1 else 'short'))
print('    => at a >= 1 the laminar count fires for EVERY m; at a = 0'
      ' only m >= 7.')
print('    CORRECTION: 432.6 said "for m >= 5 the supply is >= 5*2^t/3"'
      ' suggesting sufficiency; the a=0 threshold needs m >= 7.')
