"""multiedge410.py -- the multi-edge corridor law, the 400.2 audit, the one-edge corridor.

Colour model, per corridor edge f = (e, e+1), from 400.1 alone:
  H_f = B - {b_1,b_2} - f bipartite; O_1 - b_1 - f = two segments, chi alternates;
  400.1(ii) chi(u) = chi(u') =: c (normalize c = 0)  ==>
     chi_f(O1[p]) = p (mod 2)      for p <= e
     chi_f(O1[p]) = p + 1 (mod 2)  for p >= e+1     (A odd makes the ends agree)
  O_2 - b_2 is a path in H_f: chi_f(O2[q]) = chi_f(x) + q (mod 2).

LAW (variance): a fixed path in the middle between O_1-positions p, p' (avoiding
b_1, b_2, f - automatic) has ONE parity, which must equal
   p + p' + [p > e] + [p' > e]  (mod 2)  for EVERY corridor edge e in [s_1, s_2-1].
Non-constant over e ==> the configuration is impossible.  Same for a transition
from O_1-position p to O_2 (chi_f(O2[q]) has no e-dependence): [p > e] must be
constant over e.

CHECKS:
 (1) the 404.3 and 405.1 witness layouts have pocket attachment pairs whose
     required parity VARIES with e: both witnesses are impossible.
 (2) for every same-parity feet pair (old 400.2 model) at A in {9,11,13}, the
     interior position s_1+1 has NO legal destination: every chord, pocket
     attachment set, and transition fails the variance law.
 (3) the corrected model: adjacent feet (corridor = one edge, variance vacuous),
     mixed pin s_1+s_2+q_1+q_2 == 1 (mod 2) from P2; basic straddle room needs
     only A >= 5: the exclusion ledger's frontier is reset.
"""
import itertools

def req(p, pp, e):
    return (p + pp + (1 if p > e else 0) + (1 if pp > e else 0)) % 2

def pair_consistent(p, pp, s1, s2):
    vals = {req(p, pp, e) for e in range(s1, s2)}
    return len(vals) == 1

print('(1) witness kills')
for name, A, feet, pocket in [('404.3', 9, (2, 4), (3, 6, 7)),
                              ('405.1', 9, (1, 3), (2, 6, 7))]:
    s1, s2 = feet
    bad = [(p, pp) for p, pp in itertools.combinations(sorted(pocket), 2)
           if not pair_consistent(p, pp, s1, s2)]
    print('  witness %s  feet %s  pocket %s : inconsistent pairs %s -> %s'
          % (name, feet, pocket, bad, 'IMPOSSIBLE' if bad else 'consistent'))

print()
print('(2) interior vacancy under the old same-parity model')
for A in (9, 11, 13):
    total_cfg = 0; escapes = 0
    for s1 in range(1, A - 3):
        for s2 in range(s1 + 2, A - 2):
            if (s2 - s1) % 2: continue
            total_cfg += 1
            p = s1 + 1  # interior position
            # chords / pocket partners: any other position pp
            legal_chord = [pp for pp in range(0, A - 1)
                           if pp != p and pair_consistent(p, pp, s1, s2)]
            # pockets need >= 3 attachments incl. p: every partner must be
            # pairwise consistent with p; legal_chord empty kills both.
            # transitions: [p > e] constant over e in [s1, s2-1]?
            trans_ok = len({1 if p > e else 0 for e in range(s1, s2)}) == 1
            if legal_chord or trans_ok:
                escapes += 1
                print('  A=%d feet (%d,%d): ESCAPE chords %s trans %s'
                      % (A, s1, s2, legal_chord, trans_ok))
    print('  A=%2d  same-parity feet configs %d, interior escapes %d'
          % (A, total_cfg, escapes))

print()
print('(3) corrected form: adjacent feet, mixed pin')
for A in (5, 7, 9, 11, 13):
    ok = []
    for s1 in range(1, A - 3):
        s2 = s1 + 1
        # straddle room: t1 < s1 and t2 > s2 within 0..A-2
        if s1 >= 1 and s2 <= A - 3:
            ok.append((s1, s2))
    # mixed pin: s1 + s2 odd here (adjacent), so q1 == q2 (mod 2) is forced
    print('  A=%2d  adjacent-feet corridors with straddle room: %d %s'
          % (A, len(ok), '(ledger frontier resets to A >= 5)' if ok else ''))
