"""nc468.py -- THE (NC) TEST, RUN AT A SCALE WHERE THE ENERGY EXISTS.

466 part (E) compared the Chord Criterion against 456's counterexample
to (NC) at m = 8.  That test was VACUOUS and is recorded as such: for
girth >= 5 every pair summing to 8 needs a cycle of length <= 3, so
T_8 = 0 identically and S = 0 for every graph in the table.  Nothing
was discriminated.

The first scale with content at girth 5 is m = 16, where the splits
are 5+11, 6+10, 7+9, 8+8.  This file recomputes S at m = 16 for

  - random cubic graphs of girth >= 5 (the control),
  - 456's (NC)-breaker: a Petersen block spliced onto random blocks,
    which holds every 5-cycle inside ten vertices and drives the
    per-vertex concentration K_obs up without bound,

and asks whether K_obs and S move together.  If S stays small while
K_obs grows, the Chord Criterion is genuinely a different hypothesis
from (NC) and survives 456.  If they move together, it does not.
"""
import sys, random, itertools
from collections import Counter, deque

sys.setrecursionlimit(100000)
exec(open('tools/energy466.py').read().split(
    "ZOO = [")[0].split('"""', 2)[2])

M = 16


def profile(G, m=M, amax=6):
    """T_m by attachments, the a-profile, S, and c_z along the way."""
    cyc = all_cycles(G, m - 2)
    by_z = {}
    for c in cyc: by_z.setdefault(len(c[1]), []).append(c)
    per_a = Counter(); S = 0.0; chi1 = 0.0
    for a in range(1, amax + 1):
        z = m - 2 * a
        if z < 3: break
        cs = by_z.get(z, [])
        if not cs: continue
        R = sum(attachments(G, c, a) for c in cs)
        per_a[a] = R
        chi = R / len(cs)
        gam = 2.0 * z * len(cs) / 2 ** z
        S += chi * gam * 4.0 ** (-a) * m / z
        if a == 1: chi1 = chi
    return sum(per_a.values()), per_a, S, chi1, by_z


def Kobs(G, y):
    """456's concentration: max over v of (# y-cycles through v) over
    the mean y c_y / n."""
    n = len(G)
    cs = [c for c in all_cycles(G, y) if len(c[1]) == y]
    if not cs: return 0.0, 0
    hit = Counter()
    for c in cs:
        for v in c[0]: hit[v] += 1
    mean = y * len(cs) / n
    return max(hit.values()) / mean, len(cs)


print('=' * 72)
print('THE (NC) TEST AT m = 16  (466(E) reran; that run was vacuous)')
print('=' * 72)
print('  target 2^m/2m = %.0f' % (2.0 ** M / (2 * M)))
print()
print('  %-20s %4s %3s %8s %6s %8s %9s %9s'
      % ('graph', 'n', 'g', 'K_obs(g)', 'c_g', 'T_16', 'chi_1', 'S'))

rng = random.Random(987654321)
rows = []
for n in (20, 24, 28, 32):
    G = random_cubic(n, rng, mingirth=5)
    if G: rows.append(('random g>=5 n=%d' % n, G))

base = gp(5, 2)
for nb in (10, 14, 18, 22):
    B = random_cubic(nb, rng, mingirth=5) if nb > 10 else gp(5, 2)
    if B is None: continue
    rows.append(('Petersen+blk%d' % nb, splice(base, B)))

# a longer chain of Petersen blocks: the extreme (NC) breaker
chain = gp(5, 2)
for _ in range(2):
    chain = splice(chain, gp(5, 2))
rows.append(('Petersen chain x3', chain))

for name, G in rows:
    n = len(G); g = girth(G)
    K, cg = Kobs(G, g)
    T, per_a, S, chi1, _ = profile(G)
    print('  %-20s %4d %3d %8.2f %6d %8d %9.3f %9.5f'
          % (name, n, g, K, cg, T, chi1, S))

print()
print('  a-profiles (geometric decay in a is the claim of 466.2):')
for name, G in rows:
    T, per_a, S, chi1, _ = profile(G)
    if not T: continue
    prof = ' '.join('a=%d:%d(%.0f%%)' % (a, per_a[a], 100.0 * per_a[a] / T)
                    for a in sorted(per_a) if per_a[a])
    print('    %-20s T=%-6d %s' % (name, T, prof))

print()
print('=' * 72)
print('VERDICT')
print('=' * 72)
print('  If K_obs rises across the rows while S does not, the Chord')
print('  Criterion is a different hypothesis from (NC) and survives')
print('  456.2.  If both rise, it does not.')
