"""game414.py -- Game IV: the component census sweeps the board.

New weapon (the D_3 and pocket census, from 389.1 applied per middle
component, colours forced by the arc alternation with the f*-shift):
  colour of the attachment vertex behind O_1-position p is 1 - chi(O1[p]),
  chi(O1[p]) = (p + (1 if p > s1 else 0)) % 2.
  WINDOW CENSUS: b_3 has colour 1; beta = 1 forces every window O_1-end to
    chi = 0, so all of D_3's 1 + |WIN| deficiencies share colour 1:
        1 + |WIN| == 0 (mod 3)      (beta = 0: 1 - |WIN| == 0 (mod 3)).
  POCKET CENSUS: per pocket, (#atts with vertex-colour 1) - (#colour 0)
        == 0 (mod 3).
Soundness downgrades (Builder-favourable): D_3-jumps vetoable -> disabled;
crossing routes sterile -> impose nothing.  Kills under this weakest attack
set are unconditional.
"""
import sys, itertools, time
src = open('tools/game411.py').read().split("if __name__")[0]
exec(src)

def routes_direct(lay):
    """direct routes only: no jumps, no crossings."""
    for used, parities, trans in routes(lay):
        if trans: continue
        yield used, parities, trans

def arbitrate_weak(lay):
    A, s1, beta = lay.A, lay.s1, lay.beta
    WIN = lay.L | lay.R
    for used, parities, trans_used in routes(lay):
        if trans_used: continue                      # sterile crossings
        z = zclosure(lay, used, trans_used)
        for p in used:
            par = parities[p]
            if p in WIN:
                nb = [q for q in (p - 1, p + 1) if 0 <= q <= A - 2
                      and q not in used and not {p, q} == {s1, s1 + 1}]
                if any(q in z for q in nb):
                    if par != beta: return False
                continue
            if p in (s1, s1 + 1):
                other = s1 + 1 if p == s1 else s1
                if other not in used and other in z:
                    if par != beta: return False
                continue
            kind, atts = lay.partner[p]
            if any(a not in used and a in z for a in atts if a != p):
                if par != beta: return False
    return True

def census_ok(lay):
    A, s1, beta = lay.A, lay.s1, lay.beta
    chi = lambda p: (p + (1 if p > s1 else 0)) % 2
    WIN = lay.L | lay.R
    # window census
    wcol = [1 - chi(w) for w in WIN]                 # attachment vertex colours
    d1 = 1 + sum(1 for c in wcol if c == 1)          # b_3 is colour 1
    d0 = sum(1 for c in wcol if c == 0)
    if (d1 - d0) % 3 != 0: return False
    # pocket census
    for kind, atts in lay.part:
        if kind != 'pocket': continue
        c1 = sum(1 for a in atts if 1 - chi(a) == 1)
        c0 = len(atts) - c1
        if (c1 - c0) % 3 != 0: return False
    return True

def solve4(A, use_census=True):
    out = []
    for s1 in range(2, A - 3, 2):                    # census parity pin
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
                                if use_census and not census_ok(lay): continue
                                if arbitrate_weak(lay): out.append(lay)
    return out

if __name__ == '__main__':
    sizes = [int(x) for x in sys.argv[1:]] or [9, 11, 13, 15]
    print('weakest attack set (no jumps, sterile crossings), census imposed:')
    for A in sizes:
        t0 = time.time()
        s = solve4(A)
        print('A=%2d : %s   (%.0fs)' % (
            A, 'REFUTER WINS - EMPTY' if not s
            else 'BUILDER STANDS - %d layouts' % len(s), time.time() - t0),
            flush=True)
        for lay in s[:6]:
            print('    ', lay.desc(), flush=True)
