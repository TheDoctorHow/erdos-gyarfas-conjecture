"""game413.py -- Game III: the O_2-mirror game on the thirteen-board.

On every transition-crossing route (they exist in every Game-II layout), the
O_2-arc between the feet q_1 < q_2 is traversed; b_2's arc is Z' by
definition.  Forced certificates and their route-position parities (entry
foot at absolute parity 1, beta = 1 in every Game-II survivor):
  FEET:  q_i is interior with unused edge into the b_2-arc: certified at
         1 + m == 1 (mod 2)  ==>  m == 0: TRANSITIONS HAVE EVEN LENGTH,
         hence pi* = (A+B)/2 == 0 (mod 2), i.e. A + B == 0 (mod 4).
  INTERIOR r at arc-distance d from the entry foot: certified iff its
         structure (O_2-pocket or O_2-chord) has an attachment outside
         [q_1, q_2]; then parity 1 + m + d == 1 ==> d == 0 (mod 2).
  So odd-d interiors must sit in structures CONFINED to [q_1, q_2]; even-d
  interiors are free (certify at the right parity) or confined.
Constraints on the O_2-side layout (B odd, positions 0..B-2):
  feet same parity (mixed pin), span even >= 2; q-side third edges:
  transition-attachments only at feet-parity positions (chi pin); O_2-pockets
  >= 3 attachments (boundary law); O_2-chords span 4..B-4 (girth);
  every non-foot position must be covered.
Verdict per (B, pi*): feasible or not.
"""
import itertools

def feasible(B, pistar, verbose=False):
    n = B - 1                      # positions 0..B-2
    for par in (0, 1):             # feet parity class
        cls = [p for p in range(n) if p % 2 == par]
        for q1, q2 in itertools.combinations(cls, 2):
            span = q2 - q1
            inside = list(range(q1 + 1, q2))
            outside = [p for p in range(n) if p < q1 or p > q2]
            rest = inside + outside
            # cover rest by: transition-atts (parity par), chords, pockets
            for tset in powerset([p for p in rest if p % 2 == par]):
                rem = [p for p in rest if p not in tset]
                for part in partitions(rem, B):
                    ok = True
                    # feet certs force pi* == 0
                    if pistar != 0: ok = False
                    for kind, atts in part if ok else []:
                        confined = all(q1 <= a <= q2 for a in atts)
                        for a in atts:
                            if q1 < a < q2:
                                d = a - q1
                                if confined: continue
                                if d % 2 != 0: ok = False; break
                        if not ok: break
                    if ok:
                        if verbose:
                            print('    B=%d feet=(%d,%d) trans-atts=%s %s'
                                  % (B, q1, q2, list(tset), part))
                        return True
    return False

def powerset(xs):
    for k in range(len(xs) + 1):
        yield from itertools.combinations(xs, k)

def partitions(rem, B):
    rem = sorted(rem)
    if not rem:
        yield []
        return
    p = rem[0]; rest = rem[1:]
    for q in rest:                              # chords, girth 4..B-4
        if 4 <= q - p <= B - 4:
            for tail in partitions([x for x in rest if x != q], B):
                yield [('chord', (p, q))] + tail
    for k in range(2, len(rest) + 1):           # pockets >= 3
        for grp in itertools.combinations(rest, k):
            for tail in partitions([x for x in rest if x not in grp], B):
                yield [('pocket', (p,) + grp)] + tail

if __name__ == '__main__':
    print('O_2-side feasibility (pi* = (A+B)/2 mod 2):')
    for B in (7, 9, 11, 13):
        for pistar in (0, 1):
            f = feasible(B, pistar)
            print('  B=%2d pi*=%d : %s' % (B, pistar,
                  'feasible' if f else 'INFEASIBLE'))
    print()
    print('Thirteen-board verdict (A = 13):')
    for Bv in (7, 9, 11, 13):
        pistar = ((13 + Bv) // 2) % 2
        f = feasible(Bv, pistar)
        print('  (A,B) = (13,%2d)  pi* = %d : %s' % (Bv, pistar,
              'STANDS' if f else 'DEAD'))
