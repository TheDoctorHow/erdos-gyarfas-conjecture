"""soak472.py -- DOES THE WATER REACH? THE INTERVAL MODEL.

466.4 proved the abstract sigma-paired constant is INFINITE because
the couplings g_ij = <z_i, z_j> are free.  466.5 then showed the
object in question is the difference of two cycles meeting in a
path.  Put those together and the free coupling is exposed as a
fiction:

  A family of ears on the base cycle x has interiors disjoint
  (432.4).  Ear i has feet u_i, v_i on x and path length d_i, and its
  two arc-lifts are the cycles of length a_i + d_i and (L-a_i) + d_i
  where a_i is the chosen arc.  Two ears meet only along x, so

      Z_i cap Z_j  =  A_i cap A_j ,

  an INTERVAL OVERLAP, and the difference cycle has length

      |Z_i D Z_j|  =  l_i + l_j - 2 |A_i cap A_j| .

  THE DIMENSION COLLAPSE.  The abstract model hands the adversary
  k(k-1)/2 free couplings.  The geometry hands him 3k numbers: two
  foot positions and one path length per ear.  Quadratic against
  linear.  For k past a small threshold the coupling matrix is
  massively over-determined, and the adversary's freedom - the whole
  reason 466.4's constant was unbounded - is gone.

THE TEST.  Is the PAIRWISE constant finite in the interval model,
and is it O(2^t)?  If yes, Theorem IV follows from the pairwise
sub-problem alone, and with it (C), (D), (E) of 451.2.  Kill
conditions used, all from the archive:

  * single ear: neither arc-lift may have v_2 >= t+1 (progress from
    one ear) - island shallowness, 446.1;
  * pair: if A_i cap A_j is a single arc, the difference is a single
    cycle (457.1) and it must not have v_2 >= t+1.

Two arcs of a cycle can meet in TWO pieces; then the difference is
two disjoint cycles and there is no kill.  That escape is modelled
honestly.
"""
import sys, itertools, random
import numpy as np

P = print


def arcset(u, v, L):
    """edge positions of the arc from u forward to v (mod L)."""
    out = []
    i = u
    while i != v:
        out.append(i); i = (i + 1) % L
    return frozenset(out)


def pieces(S, L):
    """number of cyclic runs in a set of edge positions."""
    if not S: return 0
    n = 0
    for e in S:
        if (e - 1) % L not in S: n += 1
    return n if n else 1


def ears(L, dmin, dmax):
    """(u, v, a, d) with a = |arc u->v|; both arcs available via sigma."""
    out = []
    for u in range(L):
        for v in range(L):
            if u == v: continue
            a = (v - u) % L
            if a == 0: continue
            for d in range(dmin, dmax + 1):
                out.append((u, v, a, d))
    return out


def build(L, t, dmin, dmax, strict_min):
    """ear-arc choices surviving the single-ear kill, with their arcs."""
    MOD = 1 << (t + 1)
    items = []
    for (u, v, a, d) in ears(L, dmin, dmax):
        if strict_min and d <= L // 2: continue
        l1 = a + d; l2 = (L - a) + d
        if l1 % MOD == 0 or l2 % MOD == 0: continue     # single-ear progress
        items.append((u, v, a, d, l1, arcset(u, v, L)))
    return items


def compatible(A, B, L, MOD):
    """True if the pair yields NO kill."""
    (u1, v1, a1, d1, l1, S1) = A
    (u2, v2, a2, d2, l2, S2) = B
    if {u1, v1} & {u2, v2}: return None      # shared foot: excluded (432.4)
    I = S1 & S2
    if not I: return True                    # disjoint arcs: two cycles
    if pieces(I, L) != 1: return True        # two pieces: two cycles
    return (l1 + l2 - 2 * len(I)) % MOD != 0


def max_avoider(L, t, dmin, dmax, strict_min, restarts=400, seed=7):
    MOD = 1 << (t + 1)
    items = build(L, t, dmin, dmax, strict_min)
    n = len(items)
    if n == 0: return 0, 0, None
    rng = random.Random(seed)
    best = 0; bestset = None
    for _ in range(restarts):
        order = list(range(n)); rng.shuffle(order)
        cur = []
        for idx in order:
            ok = True
            for j in cur:
                c = compatible(items[idx], items[j], L, MOD)
                if c is None or c is False: ok = False; break
            if ok: cur.append(idx)
        if len(cur) > best:
            best = len(cur); bestset = [items[i][:4] for i in cur]
    return best, n, bestset


P('=' * 72)
P('THE INTERVAL MODEL: IS THE PAIRWISE CONSTANT FINITE?', flush=True)
P('=' * 72)
P('  L = 2^t m, m odd.  Ears with disjoint interiors and distinct')
P('  feet.  Adversary must dodge every single-ear and every pairwise')
P('  progress kill.  Reported: largest avoiding family found.')
P()
P('  %4s %3s %4s %8s %10s %12s %10s'
  % ('L', 't', 'm', 'ear opts', 'max avoid', 'vs 2^t', 'vs abstract'), flush=True)
for (L, t) in ((12, 2), (20, 2), (28, 2), (24, 3), (40, 3), (56, 3)):
    m = L >> t
    b, n, ws = max_avoider(L, t, 1, L - 1, False)
    P('  %4d %3d %4d %8d %10d %12s %10s'
      % (L, t, m, n, b, '%.2f x' % (b / (1 << t)), 'infinite'), flush=True)

P()
P('  Same, with 449.1 strict local minimality imposed (d > L/2 on')
P('  every ear - proved unconditionally for the residual world):')
P()
P('  %4s %3s %4s %8s %10s %12s'
  % ('L', 't', 'm', 'ear opts', 'max avoid', 'vs 2^t'), flush=True)
for (L, t) in ((12, 2), (20, 2), (28, 2), (24, 3), (40, 3), (56, 3)):
    m = L >> t
    b, n, ws = max_avoider(L, t, 1, L - 1, True)
    P('  %4d %3d %4d %8d %10d %12s'
      % (L, t, m, n, b, '%.2f x' % (b / (1 << t))), flush=True)

P()
P('=' * 72)
P('WHERE THE ADVERSARY HIDES', flush=True)
P('=' * 72)
b, n, ws = max_avoider(24, 3, 1, 23, False)
P('  L = 24, t = 3: largest family found has %d ears.' % b, flush=True)
if ws:
    P('  feet and path lengths (u, v, arc, d):', flush=True)
    for w in ws[:12]: P('    %s' % (w,), flush=True)
    As = [(w[2]) for w in ws]
    Ds = [(w[3]) for w in ws]
    P('  arcs a_i mod 8: %s' % sorted(set(a % 8 for a in As)), flush=True)
    P('  d_i  mod 8: %s' % sorted(set(d % 8 for d in Ds)), flush=True)
    P('  l_i  mod 16: %s'
      % sorted(set((w[2] + w[3]) % 16 for w in ws)), flush=True)

P()
P('=' * 72)
P('SUPPLY: HOW MANY EARS MUST ACTUALLY BE PRESENT', flush=True)
P('=' * 72)
P('  432.4 gives at least L/3 ears with disjoint interiors, and')
P('  446.2 says the full family composes for free.  Against the')
P('  measured avoider sizes:')
P()
P('  %4s %3s %10s %12s %10s' % ('L', 't', 'supply L/3', 'max avoid', 'verdict'),
  flush=True)
for (L, t) in ((12, 2), (20, 2), (28, 2), (24, 3), (40, 3), (56, 3)):
    b, n, ws = max_avoider(L, t, 1, L - 1, False)
    P('  %4d %3d %10d %12d %10s'
      % (L, t, L // 3, b, 'SPLITS' if L // 3 > b else 'holds'), flush=True)
