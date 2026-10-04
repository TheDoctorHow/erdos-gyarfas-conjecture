"""mono420.py -- target (ii) via the knife-edge: the Monochrome Board Theorem
and the mod-12 resonance.

Chain: 419.4 knife-edge ==> lambda vanishes on C_0 x C_0 (total isotropy)
==> phi(x) := lambda(x, O_1) + lambda(x, O_2) == 0 for every even cycle x.
phi = 1 on every chord-cycle (odd span by the colour law; straddling chords
via the through-b_1 completion) and on every mixed-colour pocket-cycle.
Hence: NO CHORDS anywhere on O_1 or O_2, and every middle component's
attachment set is chi-MONOCHROMATIC. The 389.1 census (per component, on
TOTALS n_1 + n_2 of O_1- and O_2-attachments) then gives n_1 + n_2 == 0
(mod 3) per component (att-vertices all one colour), |WIN| == 1 or 2 (mod 3)
by beta (b_3's extra deficiency).

BOOKS mod 3, per colour class c (pocket/carrier totals eliminated by the
component law; the census sees only totals, NOT per-O splits):
    [N_c - [w = c]|WIN|]  +  (B-1)/2  ==  0   (mod 3)
with N_0 = (A+1)/2, N_1 = (A-3)/2, w = (beta+1) % 2. Both colour classes
give the SAME condition (consistency check below):
    THE MONOCHROME RESONANCE:  A + B == 0 (mod 6) at beta = 0,
                               A + B == 4 (mod 6) at beta = 1.
With the two-sided posture A + B == 0 (mod 4):
    A + B == 0 (mod 12) at beta = 0;  A + B == 4 (mod 12) at beta = 1.

Checks: (1) phi = 1 on span419's instance-1 chord; (2) colour-count
formulas vs direct enumeration; (3) book consistency across both colour
classes, all beta; (4) execution of the recorded (35,17) beta = 0 frontier
seat; (5) beta = 1 seatability at A = 35 (seat415 machinery, beta forced);
(6) frontier sweep under floors + posture + resonance.
"""

# (1) phi on the span419 instance-1 chord (o3-o9), long completion through b_1
O1 = {('b1','o0')} | {('o%d'%i,'o%d'%(i+1)) for i in range(11)} | {('o11','b1')}
O2 = {('b2','p0')} | {('p%d'%i,'p%d'%(i+1)) for i in range(5)} | {('p5','b2')}
Zlong = {('o3','o9')} | {('o%d'%i,'o%d'%(i+1)) for i in range(3)} \
        | {('b1','o0')} | {('o11','b1')} | {('o%d'%i,'o%d'%(i+1)) for i in (9,10)}
n = lambda S: {frozenset(e) for e in S}
lam = lambda X, Y: len(n(X) & n(Y)) % 2
phi = (lam(Zlong, O1) + lam(Zlong, O2)) % 2
print('(1) chord through-b_1 completion: |Z| = %d even, lambda(Z,O_1) = %d, '
      'lambda(Z,O_2) = %d, phi = %d -> chords EXECUTED'
      % (len(Zlong), lam(Zlong, O1), lam(Zlong, O2), phi))
assert len(Zlong) % 2 == 0 and phi == 1

# (2) colour-count formulas
ok = True
for A in range(9, 120, 2):
    for s1 in range(2, A - 3, 2):
        chi = [(p + (1 if p > s1 else 0)) % 2 for p in range(A - 1)]
        ok &= chi.count(0) == (A + 1) // 2 and chi.count(1) == (A - 3) // 2
print('(2) N_0 = (A+1)/2, N_1 = (A-3)/2 for all (A odd, s_1 even): %s' % ok)
assert ok

# (3) book consistency: both colour classes give one resonance
def resonance(beta):
    out = set()
    for AB2 in range(3):                      # (A+B)/2 mod 3
        okc = []
        for c in (0, 1):
            w = (beta + 1) % 2
            WIN = 1 if w == 1 else 2          # window census mod 3
            # N_c + (B-1)/2 == (A+B)/2 + (1 if c==0 else -2) hmm -- do it
            # directly: N_0 - (A+1)/2 = 0, N_1 - (A-3)/2 = 0; book:
            # N_c - [w==c]*WIN + (B-1)/2 == 0 (mod 3), where
            # N_0 + (B-1)/2 == AB2 (mod 3), N_1 + (B-1)/2 == AB2 - 2.
            v = (AB2 if c == 0 else AB2 - 2) - (WIN if w == c else 0)
            okc.append(v % 3 == 0)
        if all(okc): out.add(AB2)
    return out
r0, r1 = resonance(0), resonance(1)
print('(3) allowed (A+B)/2 mod 3: beta=0 -> %s (A+B == 0 mod 6); '
      'beta=1 -> %s (A+B == 4 mod 6)' % (sorted(r0), sorted(r1)))
assert r0 == {0} and r1 == {2}

def alive(A, B, beta):
    if (A + B) % 4 != 0: return False                 # two-sided posture
    return (A + B) % 6 == (0 if beta == 0 else 4)     # resonance

# (4) the recorded frontier seat: (35, 17) at beta = 0
print('(4) recorded frontier seat (35,17) beta=0: alive = %s  -> EXECUTED;'
      ' (35,17) beta=1: alive = %s (pending beta=1 seatability)'
      % (alive(35, 17, 0), alive(35, 17, 1)))
assert not alive(35, 17, 0)

# (5) beta = 1 seatability at A = 35 (seat415 arithmetic, beta forced)
src = open('tools/seat415.py').read().split("if __name__")[0]
exec(src)
def seats(A, beta):
    need = 2 if beta == 1 else 1
    out = []
    for s1 in range(2, A - 3, 2):
        for kL in range(1, 10):
            if not side_ok(s1, (beta + 1) % 2, 0, kL): continue
            for kR in range(1, 10):
                if kL + kR < 3 or (kL + kR - need) % 3: continue
                if side_ok(A - 3 - s1, beta % 2, s1 + 2, kR):
                    out.append((s1, kL, kR))
    return out
s35 = seats(35, 1)
print('(5) beta=1 seats at A=35: %d found%s'
      % (len(s35), ('; first: s1=%d kL=%d kR=%d' % s35[0]) if s35 else
         '  -> (35,17) DEAD AT BOTH beta'))

# (6) frontier sweep
print('(6) survivors (floors A>=35, B>=17; posture; resonance; beta=1 only'
      ' if beta=1 seats exist at that A):')
found = []
for A in range(35, 70, 2):
    for B in range(17, A + 1, 2):
        for beta in (0, 1):
            if not alive(A, B, beta): continue
            if beta == 1 and A == 35 and not s35: continue
            found.append((A, B, beta))
for t in found[:10]:
    print('    (A,B,beta) = %s' % (t,))
print('    first survivor: %s' % (found[0],))
