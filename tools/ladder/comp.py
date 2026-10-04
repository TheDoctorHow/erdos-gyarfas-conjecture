# -*- coding: utf-8 -*-
"""515: the COMPLEMENT LAW in the closed ladder.

Chords are NOT a basis of the cycle space of a closed ladder: the host
cycle is a further generator, so a chord subset T has TWO alternating
completions, not one, and their lengths add to 2N + 2|T|.  (LQ3) is
stated for the OPEN ladder, where the completion is unique; the closed
ladder - the faithful J7 rung - therefore carries a strictly larger
target set.  This measures how much larger.
"""

import sys
from itertools import combinations
from nlutil import enum_matchings, closed_adj, girth


def census(n, step=1, cap=10 ** 9):
    N = n // 2
    lad = 0
    tsub = both = one = none = 0
    cycles = 0
    lawok = lawbad = 0
    idx = -1
    for m in enum_matchings(n):
        A = closed_adj(m, n)
        if girth(A) < 5:
            continue
        idx += 1
        if idx % step:
            continue
        if lad >= cap:
            break
        lad += 1
        ch = sorted(m)
        for r in range(0, N + 1):
            for T in combinations(ch, r):
                tsub += 1
                ends = sorted(v for e in T for v in e)
                if not ends:
                    one += 1
                    cycles += 1
                    continue
                k = len(ends)
                arcs = []
                for i in range(k):
                    a = ends[i]
                    b = ends[(i + 1) % k]
                    arcs.append((b - a) % n)
                good = 0
                lens = []
                for par in (0, 1):
                    use = [i for i in range(k) if i % 2 == par]
                    deg = {}
                    ed = []
                    for x, y in T:
                        ed.append((x, y))
                    for i in use:
                        a = ends[i]
                        b = ends[(i + 1) % k]
                        ed.append(('arc', a, b, arcs[i]))
                    adj2 = {}
                    for e in ed:
                        if e[0] == 'arc':
                            u, v = e[1], e[2]
                        else:
                            u, v = e
                        adj2.setdefault(u, []).append(v)
                        adj2.setdefault(v, []).append(u)
                    if any(len(adj2[v]) != 2 for v in adj2):
                        continue
                    s = ends[0]
                    seen, dq = set([s]), [s]
                    while dq:
                        x = dq.pop()
                        for y in adj2[x]:
                            if y not in seen:
                                seen.add(y)
                                dq.append(y)
                    if len(seen) != len(adj2):
                        continue
                    good += 1
                    lens.append(r + sum(arcs[i] for i in use))
                cycles += good
                if good == 2:
                    both += 1
                    if lens[0] + lens[1] == n + 2 * r:
                        lawok += 1
                    else:
                        lawbad += 1
                elif good == 1:
                    one += 1
                else:
                    none += 1
    print("2N=%2d ladders=%4d chord-subsets=%6d  both-completions=%6d "
          "one=%6d none=%6d | cycles=%6d  open-model targets=%6d "
          "ratio=%.3f | complement law %d ok %d bad"
          % (n, lad, tsub, both, one, none, cycles, both + one,
             float(cycles) / max(1, both + one), lawok, lawbad))
    sys.stdout.flush()


if __name__ == '__main__':
    census(12)
    census(14)
    census(16, step=8)
