# -*- coding: utf-8 -*-
"""516: k* - THE CHORD DEPTH OF THE POWER OF TWO.

Formula-free.  For every closed girth-5 ladder (host cycle 1..2N plus
a perfect matching), enumerate EVERY cycle, and for each cycle record
how many chords it uses.  Then

    k*(G) = min number of chords in a power-of-two cycle of G,
            or INF if G has none.

The archive's whole through-C apparatus lives at k* <= 3:
  k* = 1  one-chord cycle   (theta with C)
  k* = 2  two-chord cycle   (K_4 with C, when the pair crosses)
  k* = 3  three-chord cycle (K_{3,3} with C, when the triple crosses)
So the distribution of k* is exactly the coverage ledger, and every
ladder with k* >= 4 is a graph on which the layered approaches, taken
all together, deliver nothing.

Nothing here uses a weight formula, a completion, or a bridge lemma.
"""

import sys


def gen_closed_girth5(n):
    """DFS over perfect matchings of the cycle 1..n, pruning on closed
    girth >= 5.  A chord (a,b) is legal against the host cycle iff the
    circular distance min(b-a, n-(b-a)) >= 4, and against an earlier
    chord (c,d) iff the two chords are not both close on both sides."""
    mate = [0] * (n + 1)
    chords = []
    out = []

    def cdist(x, y):
        d = abs(x - y) % n
        return min(d, n - d)

    def rec():
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            out.append(list(chords))
            return
        for b in range(a + 1, n + 1):
            if mate[b]:
                continue
            if cdist(a, b) < 4:
                continue
            ok = True
            for c, d in chords:
                # 4-cycle a-b + chord + two host arcs: needs the two
                # chords' endpoints pairwise close.  General guard:
                # any cycle of length <= 4 through both chords.
                for x, y in ((c, d), (d, c)):
                    if cdist(a, x) + cdist(b, y) + 2 <= 4:
                        ok = False
                if not ok:
                    break
            if ok:
                mate[a] = b
                mate[b] = a
                chords.append((a, b))
                rec()
                chords.pop()
                mate[a] = 0
                mate[b] = 0

    rec()
    return out


def build(chords, n):
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    for a, b in chords:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def girth(adj):
    best = 99
    for s in adj:
        dist = {s: 0}
        par = {s: None}
        dq = [s]
        while dq:
            nxt = []
            for x in dq:
                for y in adj[x]:
                    if y == par[x]:
                        continue
                    if y in dist:
                        g = dist[x] + dist[y] + 1
                        if g < best:
                            best = g
                    else:
                        dist[y] = dist[x] + 1
                        par[y] = x
                        nxt.append(y)
            dq = nxt
    return best


def kstar(chords, n, P2):
    """min chords used by a power-of-two cycle; also the whole
    (length, chordcount) profile."""
    cs = set(chords)
    adj = build(chords, n)
    best = 99
    spectrum = {}
    vs = sorted(adj)
    for s in vs:
        stack = [(s, s, 1, 0, 1 << s)]
        # (current, start, length, chords used, visited bitmask)
        while stack:
            v, st, L, nc, msk = stack.pop()
            for w in adj[v]:
                e = (v, w) if v < w else (w, v)
                inc = 1 if e in cs else 0
                if w == st:
                    if L >= 3:
                        tot = L
                        c = nc + inc
                        if tot not in spectrum or spectrum[tot] > c:
                            spectrum[tot] = c
                        if tot in P2 and c < best:
                            best = c
                elif w > st and not (msk >> w & 1):
                    stack.append((w, st, L + 1, nc + inc,
                                  msk | (1 << w)))
    return best, spectrum


def sweep(n, cap=None, show=0):
    P2 = set()
    v = 4
    while v < 4 * n:
        P2.add(v)
        v *= 2
    lads = gen_closed_girth5(n)
    hist = {}
    wit = []
    for i, ch in enumerate(lads):
        if cap and i >= cap:
            break
        assert girth(build(ch, n)) >= 5
        k, _ = kstar(ch, n, P2)
        hist[k] = hist.get(k, 0) + 1
        if k >= 4:
            wit.append((k, ch))
    tot = sum(hist.values())
    print("2N=%2d  girth-5 closed ladders=%6d  k* histogram: %s"
          % (n, tot, sorted(hist.items())))
    for k, ch in wit[:show]:
        print("   k*=%d : %s" % (k, ch))
    sys.stdout.flush()
    return hist, wit


if __name__ == '__main__':
    for nn in (10, 12, 14, 16):
        sweep(nn, show=3)
