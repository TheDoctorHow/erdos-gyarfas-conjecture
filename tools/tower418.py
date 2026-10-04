"""tower418.py -- the tower's first rung, derived honestly (audit for 418).

Sources verified verbatim: 178.4/178.5 (valuation ascent; shallow = every ear
shift delta has v_2(delta) < t = v_2(|C|), i.e. delta != 0 mod 2^t); 193.4
(general t: every b_i-b_j path length avoids d == |A| - 2 mod 2^t, either arc);
195.1 (t = 2 table); 183.3 (the open family is indexed by t >= 2); 212.1 ((SC));
222.1 (mod-2^k filtration); 182.4 (modulus-8 phase boundary in three theories).

(a) GENERAL-t TABLE.  Forbidden classes per attachment pair with arc |A|:
        d == |A| - 2   and   d == -|A| - 2   (mod 2^t).
    t=2: even |A| forbids exactly ONE of the two even classes  -> (SC) pinning;
         odd  |A| forbids BOTH odd classes                     -> Parity Lemma.
    t=3: even |A| forbids <= 2 of FOUR even classes -> no pinning; odd |A|
         forbids 2 of 4 odd classes -> Parity Lemma FAILS.  (SC_8) as a
         congruence is dead.  Special feature: |A| == 2 (mod 4) forbids the
         whole 0-mod-4 sublattice -> partial mod-4 pinning survives at t = 3.
(b) IMMUNITY.  In the t=2 world (SC) allows even d == |A| (mod 4) only; then
    BOTH ear-cycle routings (2+d) + |A| and (2+d) + (|C|-|A|) are == 2 (mod 4):
    the whole ear sector dodges every power of two >= 4 by congruence alone.
(c) EXPOSURE.  Internally disjoint connecting-path pairs, same terminal pair:
    even-even (same (SC) class): sum == 0 (mod 4)   EXPOSED;
    odd-odd opposite classes:    sum == 0 (mod 4)   EXPOSED;
    odd-odd same class: sum == 2 (mod 4) immune; even-odd: odd, immune.
    Exposed theta cycles dodge {8,16,32,64,...} only by INTEGER value.
    Terminal distance >= 5 -> path length >= 5: first live rung is 16.
(d) FLOOR MEETS CEILING.  Conservative cyclomatic floor for the first seat
    (A,B) = (35,17), (beta=0, s1=8, kL=1, kR=3; mirror par 0, Q=(0,8)),
    against 222.3 isotropy ((SC) + spanning ==> d >= mu - 2) and 222.4's data
    (d >= mu - 2 never sampled past mu = 6, Petersen).
"""

def table(t):
    M = 1 << t
    rows = {}
    for r in range(M):
        rows[r] = sorted({(r - 2) % M, (-r - 2) % M})
    return rows

print('=== (a) general-t forbidden tables ===')
for t in (2, 3, 4):
    M = 1 << t
    rows = table(t)
    ev = [c for c in range(M) if c % 2 == 0]
    od = [c for c in range(M) if c % 2 == 1]
    pin_even = all(len([c for c in ev if c not in rows[r]]) == 1
                   for r in range(0, M, 2))
    parity_lemma = all(all(c in rows[r] for c in od)
                       for r in range(1, M, 2))
    print(' t=%d (mod %d):' % (t, M))
    for r in range(M):
        print('   |A|=%2d  forbidden %s' % (r, rows[r]))
    print('   single-even-class pinning at even |A|: %s' % pin_even)
    print('   Parity Lemma (odd |A| kills all odd classes): %s' % parity_lemma)
    if t == 3:
        for r in (2, 6):
            allow_ev = [c for c in ev if c not in rows[r]]
            assert allow_ev == [2, 6], allow_ev
        print('   t=3 special: |A| == 2 (mod 4) confines evens to 2 mod 4: True')

print()
print('=== (b) immunity of the ear sector (t=2 world) ===')
fails = checks = 0
for Cl in range(12, 104, 4):                     # |C| == 0 mod 4
    for Aarc in range(2, Cl - 1, 2):             # arcs even (same-parity world)
        c = Aarc % 4                             # allowed even class = |A| mod 4
        for d in range(c, 101, 4):
            for E in ((2 + d) + Aarc, (2 + d) + (Cl - Aarc)):
                checks += 1
                if E % 4 != 2:
                    fails += 1
print(' ear-cycles == 2 mod 4: %d checks, %d failures' % (checks, fails))

print()
print('=== (c) exposure of theta pairs ===')
def first_hits(sums, tag):
    pows = [p for p in (8, 16, 32, 64, 128) if p in sums]
    print(' %s: all sums == 0 mod 4: %s; powers of 2 realizable: %s'
          % (tag, all(s % 4 == 0 for s in sums), pows))
    for p in pows[:2]:
        wit = sorted(sums[p])[:4]
        print('   %d = %s' % (p, ' = '.join('%d+%d' % w for w in wit)))
sums_ee = {}
for c in (0, 2):
    ds = [d for d in range(6, 130, 2) if d % 4 == c]
    for i, d in enumerate(ds):
        for dd in ds[i:]:
            sums_ee.setdefault(d + dd, set()).add((d, dd))
first_hits(sums_ee, 'even-even same class (forced by (SC))')
sums_oo = {}
for d in range(5, 130, 2):
    for dd in range(5, 130, 2):
        if d % 4 != dd % 4 and d <= dd:
            sums_oo.setdefault(d + dd, set()).add((d, dd))
first_hits(sums_oo, 'odd-odd opposite classes')
im = [(d + dd) % 4 for d in range(5, 60, 2) for dd in range(5, 60, 2)
      if d % 4 == dd % 4]
print(' odd-odd same class: sums == 2 mod 4 always: %s (immune)'
      % all(x == 2 for x in im))

print()
print('=== (d) frontier floor vs isotropy ceiling ===')
A, B = 35, 17
WIN, Q = 4, 2                                    # first seats: kL+kR = 4; |Q| = 2
items = [
    ('O1 structures  (%d covered positions, >= 1 cycle per 2 atts)'
     % ((A - 1) - 2 - WIN), ((A - 1) - 2 - WIN) // 2),
    ('D3 component   (attaches at %d windows)' % WIN, WIN - 1),
    ('corridor cycle (f*)', 1),
    ('O2 structures  (%d covered positions)' % ((B - 1) - Q), ((B - 1) - Q) // 2),
    ('O2 transitions (%d feet)' % Q, Q - 1),
    ('b_1, b_2 restored (degree >= 2 each)', 2),
]
mu = 0
for name, v in items:
    mu += v
    print('   %-58s %2d' % (name, v))
print('   conservative cyclomatic FLOOR ESTIMATE: mu(B) >= %d' % mu)
print('   isotropy demand under (SC)+spanning: d >= mu - 2 >= %d' % (mu - 2))
print('   largest sampled mu with d >= mu-2 (222.4, ~2600 graphs): 6 (Petersen)')
lo_pair = 5 + 7                                  # smallest exposed odd-odd sum
print('   live rungs at frontier (pair sums from %d up, reach ~A+B=%d+): '
      '16, 32%s' % (lo_pair, A + B, ', 64 once |B| ~ 64 -- the escape direction'))
