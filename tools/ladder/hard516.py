# -*- coding: utf-8 -*-
"""516: the through-C template ledger, measured in the hard case.

Open girth-legal ladder on [1, 2N]: perfect matching, every span >= 4,
every non-disjoint pair at L1 >= 3.  HARD CASE (501.0): additionally
no span in {7, 15, 31} and no non-disjoint pair with L1 in {6, 14, 30}.

Templates through C in the closed ladder (515.3 + 514.6):
  C + chord           = theta:  2N, s+1, 2N+1-s
  C + crossing pair   = K_4:    + 2+L1, 2N-L1
  C + crossing triple = K_{3,3}: + XXX, 2N+6-XXX
In the hard case every theta and K_4 value is a non-power by fiat
(2N never a power at the sizes run), so the whole template machinery
through C delivers exactly XXX and 2N+6-XXX per crossing triple.

MEASURED here, exhaustively at 2N = 14..22:
  col A: hard-case ladder count            (cross-check 501.0)
  col B: ladders where no crossing triple has XXX a power of two
                                           (cross-check 508.2: 33, 24)
  col C: ladders where no crossing triple has XXX or 2N+6-XXX a power
         of two  -- THE NEW NUMBER: the residue of the full
         through-C template ledger in the hard case
  col D: of col C, how many are also hard in the CLOSED reading
         (no closed one-chord value 2N+1-s and no closed two-chord
         value -- crossing complement 2N-L1, disjoint-pair wrap --
         a power of two)
  col E: of col D, how many have closed girth >= 5
For every col-C ladder, verify (TRI): some three-chord open value
3 + sum of block gaps is a power of two (501.2 criterion: no block
{q1,q2},{q3,q4},{q5,q6} is itself a chord).
"""

import sys

P2 = set()
v = 4
while v < 300:
    P2.add(v)
    v *= 2

BADSPAN = {7, 15, 31}
BADL1 = {6, 14, 30}


def gen_hard(n, out):
    """DFS over matchings of [1,n]; open girth-legal + hard case."""
    mate = [0] * (n + 1)
    chords = []

    def rec():
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            out.append(list(chords))
            return
        for b in range(a + 4, n + 1):
            if mate[b]:
                continue
            s = b - a
            if s in BADSPAN:
                continue
            ok = True
            for c, d in chords:
                if d > a:  # non-disjoint (c < a always)
                    l1 = (a - c) + abs(b - d)
                    if l1 < 3 or l1 in BADL1:
                        ok = False
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


def closed_girth5(chords, n):
    """exact girth >= 5: for every edge uv, the shortest u-v path
    avoiding uv must have length >= 4 (BFS)."""
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    for a, b in chords:
        adj[a].append(b)
        adj[b].append(a)
    edges = set()
    for v in adj:
        for w in adj[v]:
            edges.add((v, w) if v < w else (w, v))
    for u, v in edges:
        dist = {u: 0}
        dq = [u]
        found = False
        while dq and not found:
            nxt = []
            for x in dq:
                if dist[x] >= 3:
                    break
                for y in adj[x]:
                    if x == u and y == v:
                        continue
                    if y not in dist:
                        dist[y] = dist[x] + 1
                        if y == v:
                            found = True
                        nxt.append(y)
            dq = nxt
        if found and dist[v] <= 3:
            return False
    return True


def tri_ok(chords, n):
    """(TRI): some legal 3-chord open completion has power-of-2 weight."""
    cs = set(chords)
    from itertools import combinations
    for t in combinations(chords, 3):
        q = sorted(v for e in t for v in e)
        if (q[0], q[1]) in cs or (q[2], q[3]) in cs or (q[4], q[5]) in cs:
            continue
        w = 3 + (q[1] - q[0]) + (q[3] - q[2]) + (q[5] - q[4])
        if w in P2:
            return True
    return False


def sweep(n):
    ladders = []
    gen_hard(n, ladders)
    nb = nc = nd = ne = 0
    tri_fail = []
    residues = []
    for ch in ladders:
        ch = sorted(ch)
        k = len(ch)
        hitx = hitc = False
        for i in range(k):
            a1, b1 = ch[i]
            if hitc and hitx:
                break
            for j in range(i + 1, k):
                a2, b2 = ch[j]
                if a2 >= b1:
                    break
                for l in range(j + 1, k):
                    a3, b3 = ch[l]
                    if a3 >= b1:
                        break
                    if not (b1 < b2 < b3 and a3 < b1):
                        continue
                    xxx = 3 + (a2 - a1) + (b1 - a3) + (b3 - b2)
                    if xxx in P2:
                        hitx = True
                        hitc = True
                    if (n + 6 - xxx) in P2:
                        hitc = True
        if not hitx:
            nb += 1
        if not hitc:
            nc += 1
            if not tri_ok(ch, n):
                tri_fail.append(ch)
            # col D: closed-legal (no span above 2N-4; closed girth 5)
            if any(b - a > n - 4 for a, b in ch):
                continue
            if not closed_girth5(ch, n):
                continue
            nd += 1
            # col E: the closed one/two-chord layer misses too.
            # values: 2N+1-s per chord; 2N+2-L1 per CROSSING pair
            # (both completions exist only there); the single wrap
            # value per line-disjoint pair; nested pairs add nothing
            # (their one completion 2+L1 is hard-banned).
            chit = False
            for a, b in ch:
                if (n + 1 - (b - a)) in P2:
                    chit = True
            for i in range(k):
                a1, b1 = ch[i]
                for j in range(i + 1, k):
                    a2, b2 = ch[j]
                    if a2 < b1 and b2 > b1:  # crossing
                        l1 = (a2 - a1) + (b2 - b1)
                        if (n + 2 - l1) in P2:
                            chit = True
                    elif a2 > b1:  # line-disjoint: one wrap cycle
                        w = 2 + (a2 - b1) + (n - b2 + a1)
                        if w in P2:
                            chit = True
            if not chit:
                ne += 1
                residues.append(ch)
    print("2N=%2d  A hard=%6d  B openXXX-fail=%5d  C triple-fail=%5d  "
          "D closed-legal=%5d  E RESIDUE=%4d  TRI-fail=%d"
          % (n, len(ladders), nb, nc, nd, ne, len(tri_fail)))
    for ch in tri_fail:
        print("   TRI FAILURE:", ch)
    for ch in residues[:6]:
        print("   residue:", ch)
    sys.stdout.flush()
    return len(ladders), nb, nc, nd, ne, residues


if __name__ == '__main__':
    for n in (14, 16, 18, 20, 22):
        sweep(n)
