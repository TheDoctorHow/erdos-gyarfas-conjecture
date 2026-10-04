"""game412.py -- Game II: the census-pinned endgame.

Game I (tools/game411.py) omitted a standing law: the census forcing 400.1(iv)
read in corridor coordinates.  chi(u) = chi(u') = 0 and chi alternates to g at
position s_1 (chi(g) = s_1 mod 2), and the mod-3 census forces the quadruple
u, u', g, h into ONE colour class (H is connected: {b_1,b_2} is no cut by the
2-Cut Lemma, and f* is no bridge by the window straddle):
        CENSUS PARITY PIN:   s_1 == 0  (mod 2).
Every Game-I survivor has s_1 odd: all are executed.  Game II re-arbitrates the
whole board with the pin imposed, at A = 9..17.
"""
import sys
src = open('tools/game411.py').read()
src = src.split("if __name__")[0]
exec(src)

def solve2(A):
    survivors = []
    for s1 in range(2, A - 3, 2):          # census parity pin: s_1 even
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
    return survivors

if __name__ == '__main__':
    import time
    sizes = [int(x) for x in sys.argv[1:]] or [9, 11, 13, 15]
    for A in sizes:
        t0 = time.time()
        s = solve2(A)
        print('A=%2d : %s   (%.0fs)' % (
            A, 'REFUTER WINS (census-pinned board empty)' if not s
            else 'BUILDER STANDS - %d layouts' % len(s), time.time() - t0),
            flush=True)
        for lay in s[:8]:
            print('    ', lay.desc(), flush=True)
