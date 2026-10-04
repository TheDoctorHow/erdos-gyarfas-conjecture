"""game411.py -- the certificate game on the corrected canonical form, arbitrated.

The {1,3}-monochromy endgame (sealed template world, one-edge corridor f* =
(s_1, s_1+1), Law R: a {1,3}-connecting path is even iff it avoids f*).
Builder chooses beta, feet, window set WIN (L left of f*, R right), and a
partition of remaining positions into chords (4 <= span <= A-4) and pockets
(>= 3 attachments).  Refuter plays every route skeleton:
  b_1 enters at position 0 or A-2, walks arc segments, may jump through D_3
  at window attachments (each used once, <= 2 jumps), may cross once through
  the transitions (s_1 <-> s_1+1 via O_2), and exits into D_3 at a window
  attachment (to b_3).  All f*-avoiding routes are EVEN (Law R), so every one
  is in play.
Certification (Refuter-conservative, all forced):  an interior position p is
certified iff its third-edge partner lies in the forced Z'-closure of b_2:
seeded by O_2 (the b_2-arc is never used), grown through unused transitions to
their feet, along unused-arc adjacency (including f* if both feet unused), and
through structures with an attachment already in the closure.  Window third
edges (into D_3) and partners on the route are safe for the Builder.
Builder survives a layout iff EVERY certified position on EVERY route has one
route-position parity == beta.  Output per A: surviving layouts.
"""
import itertools, sys

def solve(A, verbose=False):
    survivors = []
    for s1 in range(1, A - 3):
        feet = (s1, s1 + 1)
        avail = [p for p in range(A - 1) if p not in feet]
        low = [p for p in avail if p < s1]
        high = [p for p in avail if p > s1 + 1]
        for beta in (0, 1):
            Lc = [p for p in low if p % 2 == (beta + 1) % 2]
            Rc = [p for p in high if p % 2 == beta % 2]
            for ls in range(1, len(Lc) + 1):
                for L in itertools.combinations(Lc, ls):
                    for rs in range(1, len(Rc) + 1):
                        for R in itertools.combinations(Rc, rs):
                            if len(L) + len(R) < 3: continue
                            WIN = set(L) | set(R)
                            rem = [p for p in avail if p not in WIN]
                            for part in partitions(rem, A):
                                lay = Layout(A, s1, beta, set(L), set(R), part)
                                if arbitrate(lay):
                                    survivors.append(lay)
                                    if verbose:
                                        print('   SURVIVES', lay.desc())
    return survivors

def partitions(rem, A):
    """all ways to cover rem by chords (pairs, 4<=span<=A-4) and pockets (>=3)."""
    rem = sorted(rem)
    if not rem:
        yield []
        return
    p = rem[0]
    rest = rem[1:]
    # chord partners
    for q in rest:
        if 4 <= q - p <= A - 4:
            for tail in partitions([x for x in rest if x != q], A):
                yield [('chord', (p, q))] + tail
    # pockets containing p: choose >= 2 partners
    for k in range(2, len(rest) + 1):
        for grp in itertools.combinations(rest, k):
            for tail in partitions([x for x in rest if x not in grp], A):
                yield [('pocket', (p,) + grp)] + tail

class Layout:
    def __init__(self, A, s1, beta, L, R, part):
        self.A, self.s1, self.beta, self.L, self.R, self.part = A, s1, beta, L, R, part
        self.partner = {}
        for kind, atts in part:
            for a in atts:
                self.partner[a] = (kind, atts)
    def desc(self):
        return 'A=%d feet=(%d,%d) beta=%d L=%s R=%s %s' % (
            self.A, self.s1, self.s1 + 1, self.beta,
            sorted(self.L), sorted(self.R), self.part)

def routes(lay):
    """yield (usedset, posparity: dict pos -> route-position parity, trans_used)"""
    A, s1 = lay.A, lay.s1
    WIN = lay.L | lay.R
    def chi(p):  # colour of O1[p], f* between s1, s1+1
        return (p + (1 if p > s1 else 0)) % 2
    def walk(pos, par, used, parities, jumps, trans_used):
        # at O_1 position pos with route-parity par (pos already recorded)
        # option 1: exit here if window attachment -> complete route
        if pos in WIN:
            yield used, dict(parities), trans_used
        # option 2: continue along arc both directions
        for step in (-1, 1):
            nxt = pos + step
            if nxt < 0 or nxt > A - 2: continue
            if {pos, nxt} == {s1, s1 + 1}: continue      # f* forbidden
            if nxt in used: continue
            u2 = used | {nxt}; parities[nxt] = (par + 1) % 2
            yield from walk(nxt, (par + 1) % 2, u2, parities, jumps, trans_used)
            del parities[nxt]
        # option 3: D_3 jump from a window attachment
        if pos in WIN and jumps < 2:
            for t in WIN:
                if t in used: continue
                lam = (chi(pos) + chi(t)) % 2   # forced D_3-path parity
                p2 = (par + 1 + lam + 1) % 2
                u2 = used | {t}; parities[t] = p2
                yield from walk(t, p2, u2, parities, jumps + 1, trans_used)
                del parities[t]
        # option 4: transition crossing (once), from a foot to the other foot
        if pos in (s1, s1 + 1) and not trans_used:
            other = s1 + 1 if pos == s1 else s1
            if other not in used:
                p2 = (par + 0) % 2   # m1 + |q2-q1| + m2 == 0 (mod 2)
                u2 = used | {other}; parities[other] = p2
                yield from walk(other, p2, u2, parities, 2 if False else 0 + (jumps), True)
                del parities[other]
    for start in (0, A - 2):
        used = {start}; parities = {start: 1}
        yield from walk(start, 1, used, parities, 0, False)

def zclosure(lay, used, trans_used):
    """forced Z'(route) closure among unused positions; O_2-side always seeded."""
    A, s1 = lay.A, lay.s1
    unused = set(range(A - 1)) - used
    z = set()
    seed = []
    if not trans_used:
        for f in (s1, s1 + 1):
            if f in unused: seed.append(f)   # foot -> its transition -> O_2
    frontier = list(seed)
    z |= set(seed)
    changed = True
    while changed:
        changed = False
        for p in list(unused - z):
            # arc adjacency (f* usable: s1~s1+1 edge exists)
            for q in (p - 1, p + 1):
                if q in z and 0 <= q <= A - 2:
                    z.add(p); changed = True; break
            if p in z: continue
            # structure link
            if p in lay.partner:
                kind, atts = lay.partner[p]
                if any(a in z for a in atts if a != p):
                    z.add(p); changed = True
    return z

def arbitrate(lay):
    A, s1, beta = lay.A, lay.s1, lay.beta
    WIN = lay.L | lay.R
    for used, parities, trans_used in routes(lay):
        z = zclosure(lay, used, trans_used)
        if not trans_used:
            pass  # O_2 reachable only via unused feet (in closure seed)
        for p in used:
            par = parities[p]
            # third-edge target of p:
            if p in WIN:
                # exit/window: D_3 target killable by Builder EXCEPT the final
                # exit's unused O_1 edge is NOT this; window att interior via
                # both O_1 edges has D_3 third edge: skip (Builder controls D_3)
                # BUT if p was the route's exit, its unused O_1-edge neighbour:
                # exits recorded as window: conservative skip here; the unused
                # O_1-edge case is covered because exit vertices in `used` with
                # one O_1-neighbour unused certify via arc adjacency below.
                nb = [q for q in (p - 1, p + 1) if 0 <= q <= A - 2
                      and q not in used and not {p, q} == {s1, s1 + 1}]
                if any(q in z for q in nb):
                    if par != beta: return False
                continue
            if p in (s1, s1 + 1):
                other = s1 + 1 if p == s1 else s1
                tgt = None
                if other not in used and other in z:
                    if par != beta: return False
                continue
            kind, atts = lay.partner[p]
            others = [a for a in atts if a != p]
            if any(a not in used and a in z for a in others):
                if par != beta: return False
    return True

if __name__ == '__main__':
    for A in (5, 7, 9, 11):
        s = solve(A)
        print('A=%2d : %s' % (A, 'REFUTER WINS (no layout survives)' if not s
                              else 'BUILDER STANDS - %d layouts' % len(s)))
        seen = set()
        for lay in s[:6]:
            d = lay.desc()
            if d not in seen: seen.add(d); print('    ', d)
