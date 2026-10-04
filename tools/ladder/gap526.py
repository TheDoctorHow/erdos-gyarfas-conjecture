# -*- coding: utf-8 -*-
"""526: is (ATOM) true?  (ATOM) says the spectrum meets [x, 2x] for
every x with g <= 2x <= c.  Equivalently, listing the spectrum as
l_1 < ... < l_r, it says NO CONSECUTIVE PAIR HAS l_{i+1} >= 2 l_i + 2.

So (ATOM) is a NO-RATIO-2-GAP statement, and it is refuted by any
graph with an isolated short cycle.  This measures the largest gap
l_{i+1} - 2 l_i over random girth-5 closed ladders, looking for one
that reaches +2.
"""

import sys
import random
from itertools import combinations


def girth_legal_closed(ch, n):
    def cd(x, y):
        d = abs(x - y) % n
        return min(d, n - d)
    for a, b in ch:
        if cd(a, b) < 4:
            return False
    for (a, b), (c, d) in combinations(ch, 2):
        for x, y in ((c, d), (d, c)):
            if cd(a, x) + cd(b, y) + 2 <= 4:
                return False
    return True


def spectrum(ch, n, cap):
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    for a, b in ch:
        adj[a].append(b)
        adj[b].append(a)
    got = set()
    vs = sorted(adj)
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        st = [(s, (s,), 1 << s)]
        while st:
            v, path, msk = st.pop()
            for w in adj[v]:
                if w == s:
                    if len(path) >= 3:
                        got.add(len(path))
                elif (rank[w] > rank[s] and not (msk >> w & 1)
                      and len(path) < cap):
                    st.append((w, path + (w,), msk | (1 << w)))
    return sorted(got)


def rand_ladder(n, rnd, tries=3000):
    for _ in range(tries):
        pts = list(range(1, n + 1))
        rnd.shuffle(pts)
        ch = [tuple(sorted((pts[2 * i], pts[2 * i + 1])))
              for i in range(n // 2)]
        if girth_legal_closed(ch, n):
            return ch
    return None


def worst_gap(n, samples, cap):
    rnd = random.Random(4242 + n)
    best = None
    for _ in range(samples):
        ch = rand_ladder(n, rnd)
        if ch is None:
            continue
        sp = spectrum(ch, n, cap)
        if len(sp) < 2:
            continue
        for i in range(len(sp) - 1):
            d = sp[i + 1] - 2 * sp[i]
            if best is None or d > best[0]:
                best = (d, sp[i], sp[i + 1], list(ch), sp[:8])
    print("2N=%2d  worst  l_{i+1} - 2 l_i = %+d   at (%d -> %d)   "
          "%s" % (n, best[0], best[1], best[2],
                  "*** (ATOM) FAILS" if best[0] >= 2 else "no gap"))
    print("       spectrum head:", best[4])
    sys.stdout.flush()
    return best


if __name__ == '__main__':
    for n, s, c in ((20, 300, 21), (30, 150, 20), (40, 60, 18)):
        worst_gap(n, s, c)
