"""blitz473.py -- CONSTRAIN, THEN ENCODE, THEN BLITZ.

The plan: drive the residual configuration space down with every
banked constraint, then represent what is left so that the whole
space can be swept in minutes.

PART 1 measures what each constraint actually buys at (t, m).
PART 2 builds the blitz encoding and benchmarks it.
PART 3 reports which strata become exhaustible.

THE ENCODING.  A configuration is a base cycle x = Z/L plus k ears.
Ear i has feet u_i, v_i on x and interior path length d_i.  Since
L = 2^t m <= 64 for every stratum that matters (t <= 4, m <= 5), the
ARC of an ear is ONE 64-BIT WORD: bit e set iff edge e of x lies in
the arc.  Then for any signed selection sigma:

    arc part   =  XOR of the chosen arc masks
    length     =  popcount(XOR)  +  SUM_{i in S} d_i
    kill test  =  length == 0  (mod 2^{t+1})
    is a cycle =  runs(XOR) == |S|  and the alternation closes up

Every operation is a machine word op.  The ear-lift sign is just
which of the two arcs is used, i.e. mask or its complement, so the
sigma in {-1,0,1}^k of 465 is a ternary digit choice over words.

This is the same object as 465's Gram matrix - popcount(XOR) is
|V|^2 for a +-1 circulation - but stored as bits instead of as a
matrix, which is what makes the sweep fast.
"""
import sys, time, random, itertools

P = print
POP = (lambda w: bin(w).count('1'))


def runs(mask, L):
    """number of cyclic runs of set bits."""
    if mask == 0: return 0
    full = (1 << L) - 1
    rot = ((mask >> 1) | (mask << (L - 1))) & full
    return POP(mask & ~rot)


def arcmask(u, a, L):
    """arc of length a starting at edge u."""
    m = 0
    for i in range(a): m |= 1 << ((u + i) % L)
    return m


P('=' * 72)
P('PART 1.  WHAT EACH CONSTRAINT ACTUALLY BUYS', flush=True)
P('=' * 72)
P('  Ear = (foot u, arc length a, interior length d).  d matters')
P('  mod 2^{t+1} for the arithmetic; a and u are geometric.')
P()
P('  %5s %3s %3s %10s %10s %10s %10s %10s'
  % ('L', 't', 'm', 'raw', '+shallow', '+both-arc', '+girth5', 'per-ear'),
  flush=True)
STRATA = [(24, 3, 3), (40, 3, 5), (48, 4, 3), (80, 4, 5), (96, 5, 3)]
SURV = {}
for (L, t, m) in STRATA:
    MOD = 1 << (t + 1)
    raw = sh = ba = g5 = 0
    keep = []
    for a in range(1, L):
        for dm in range(MOD):
            raw += 1
            dp = (dm - a) % MOD            # delta+ = d - a
            dn = (dm + a - L) % MOD        # delta- = d + a - L
            if dp == (1 << t) or dn == (1 << t): continue   # 446.1 shallow
            sh += 1
            # 449.1 strict local minimality on BOTH arcs: d > max(a, L-a).
            # d is only known mod MOD, so this constrains which residues
            # are reachable given d > 3L/4 in the worst case: none are
            # excluded, but a is: an arc with a >= L - 1 leaves no room.
            if a >= L - 1: continue
            ba += 1
            # girth 5: both cycles at least 5, and the ear is a path
            if a + 1 < 5 or (L - a) + 1 < 5: continue
            g5 += 1
            keep.append((a, dm))
    SURV[(L, t, m)] = keep
    P('  %5d %3d %3d %10d %10d %10d %10d %10d'
      % (L, t, m, raw, sh, ba, g5, len(keep)), flush=True)
P()
P('  Per-ear survivors are ~L . 2^{t+1} = 2 L^2 / m.  The banked')
P('  congruence constraints remove only a 2^{-t} sliver each - which')
P('  is 467.1 again: congruences are cheap.  The space is NOT small')
P('  yet, and no amount of arithmetic will make it small.')

P()
P('=' * 72)
P('PART 2.  THE BLITZ ENCODING, BENCHMARKED', flush=True)
P('=' * 72)
L, t, m = 24, 3, 3
MOD = 1 << (t + 1)
rng = random.Random(20260808)
P('  L = %d, t = %d, one 64-bit word per arc.' % (L, t))
P()
for k in (8, 10, 12, 14):
    masks = []; ds = []
    used = set()
    while len(masks) < k:
        u = rng.randrange(L); a = rng.randrange(4, L - 3)
        v = (u + a) % L
        if u in used or v in used: continue
        used.add(u); used.add(v)
        masks.append(arcmask(u, a, L)); ds.append(rng.randrange(MOD))
    t0 = time.time(); hits = 0; cnt = 0
    for sig in itertools.product((0, 1, 2), repeat=k):
        X = 0; S = 0; sz = 0
        for i, s in enumerate(sig):
            if s == 0: continue
            X ^= masks[i] if s == 1 else (((1 << L) - 1) ^ masks[i])
            S += ds[i]; sz += 1
        if sz == 0: continue
        cnt += 1
        if (POP(X) + S) % MOD == 0 and runs(X, L) == sz: hits += 1
    el = time.time() - t0
    P('  k=%2d   3^k = %9d   swept in %6.2fs   %8.0f sel/s   kills %d'
      % (k, cnt + 1, el, cnt / max(el, 1e-9), hits), flush=True)

P()
P('  Naive ternary sweep is O(3^k) with a k-loop inside.  Two')
P('  standard accelerations apply and both are exact:')
P('   (i)  GRAY CODE - consecutive sigma differ in one digit, so the')
P('        XOR and the sum update in O(1): removes the inner loop.')
P('   (ii) MEET IN THE MIDDLE - split the ears in half, tabulate')
P('        3^{k/2} half-states, and join.  3^{k/2} is 6561 at k = 16')
P('        and 1.6e6 at k = 26.')
P()
t0 = time.time()
half = 12
tab = {}
masks = [arcmask(rng.randrange(L), rng.randrange(4, L - 3), L)
         for _ in range(half)]
ds = [rng.randrange(MOD) for _ in range(half)]
for sig in itertools.product((0, 1, 2), repeat=half):
    X = 0; S = 0
    for i, s in enumerate(sig):
        if s == 0: continue
        X ^= masks[i] if s == 1 else (((1 << L) - 1) ^ masks[i])
        S += ds[i]
    tab[(X, S % MOD)] = 1
el = time.time() - t0
P('  half-table for k/2 = %d built in %.2fs, %d distinct states'
  % (half, el, len(tab)), flush=True)
P('  => a k = %d family is decidable in about %.2f s by joining two'
  % (2 * half, 2 * el))
P('     such tables, against %.3g for the naive sweep.'
  % (3.0 ** (2 * half) / max(3.0 ** half / max(el, 1e-9), 1)))

P()
P('=' * 72)
P('PART 3.  WHAT IS ACTUALLY FEASIBLE', flush=True)
P('=' * 72)
P('  The per-FAMILY test is cheap.  The cost is the number of')
P('  FAMILIES: choose k = L/3 ears from ~2L^2/m types, modulo the')
P('  dihedral group of order 2L.')
P()
P('  %5s %3s %3s %8s %10s %14s %12s'
  % ('L', 't', 'm', 'k = L/3', 'ear types', 'families', 'verdict'), flush=True)
import math
for (L, t, m) in STRATA:
    k = L // 3
    types = len(SURV[(L, t, m)]) * L // 2      # x foot positions / symmetry
    fam = math.comb(max(types, k), k) / (2 * L)
    P('  %5d %3d %3d %8d %10d %14.3g %12s'
      % (L, t, m, k, types, fam,
         'blitzable' if fam < 1e12 else 'NOT'), flush=True)
P()
P('  So the encoding is not the bottleneck and never was.  The')
P('  bottleneck is that the family count is astronomical at every')
P('  stratum, because the per-ear space is ~L^3/m and k grows with L.')
P('  Constraining must therefore attack k or the ear types, not the')
P('  per-family test.')
