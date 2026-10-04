# -*- coding: utf-8 -*-
"""518: (E_3) at girth 6 and 7 - the two bands-3 girths never searched.

(E_3): every cubic graph of girth 5, 6 or 7 has an 8-cycle.  516.4 and
517.2 searched girth 5.  Here girth >= 6 and girth >= 7, over
Hamiltonian cubic graphs (= closed ladders), with two monotone prunes:
no cycle of length < gmin, and no cycle of length 8.  Both survive
edge addition, so absence is proved on the class.
"""

import sys


def search(n, gmin, want=2):
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    mate = [0] * (n + 1)
    chords = []
    hits = []
    nodes = [0]

    def cdist(x, y):
        d = abs(x - y) % n
        return min(d, n - d)

    def path_exists(a, b, L, ban):
        stack = [(a, 0, 1 << a)]
        while stack:
            v, d, msk = stack.pop()
            for w in adj[v]:
                if frozenset((v, w)) in ban:
                    continue
                if w == b:
                    if d + 1 == L:
                        return True
                    continue
                if d + 1 < L and not (msk >> w & 1):
                    stack.append((w, d + 1, msk | (1 << w)))
        return False

    def ok(a, b):
        # a path of length L plus the chord closes a cycle of length
        # L+1, so girth >= gmin bans L+1 <= gmin-1, i.e. L <= gmin-2.
        # The old range(2, gmin) banned L+1 = gmin as well and so
        # enforced girth >= gmin+1.  CORRECTED.
        ban = set([frozenset((a, b))])
        for L in range(2, gmin - 1):
            if path_exists(a, b, L, ban):
                return False
        return not path_exists(a, b, 7, ban)

    def rec():
        nodes[0] += 1
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            hits.append(list(chords))
            return
        for b in range(a + 1, n + 1):
            if mate[b] or cdist(a, b) < gmin - 2:
                continue
            adj[a].append(b)
            adj[b].append(a)
            if ok(a, b):
                mate[a] = b
                mate[b] = a
                chords.append((a, b))
                rec()
                chords.pop()
                mate[a] = 0
                mate[b] = 0
            adj[a].remove(b)
            adj[b].remove(a)
            if len(hits) >= want:
                return

    rec()
    print("girth>=%d  n=%2d  Hamiltonian cubic, no 8-cycle: %s%d  "
          "(nodes %d)" % (gmin, n, ">=" if len(hits) >= want else "",
                          len(hits), nodes[0]))
    for ch in hits:
        print("    ", ch)
    sys.stdout.flush()
    return hits


if __name__ == '__main__':
    import sys as S
    gm = int(S.argv[1])
    lo, hi = int(S.argv[2]), int(S.argv[3])
    for nn in range(lo, hi + 1, 2):
        search(nn, gm)
