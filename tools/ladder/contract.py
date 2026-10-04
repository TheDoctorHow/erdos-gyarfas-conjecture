"""Chord contraction on a closed ladder.

Delete a chord {u,v} and suppress its two endpoints: u's host neighbours
u-1, u+1 become adjacent, likewise v-1, v+1.  The host cycle survives as
a cycle of length 2N-2 and the remaining N-1 chords still form a perfect
matching, so the result is again a closed ladder.
"""
import sys, random
from collections import deque

def contract(mate, n, u):
    """contract the chord at u.  Returns the new mate array on n-2 vertices."""
    v = mate[u]
    dead = {u, v}
    keep = [x for x in range(n) if x not in dead]
    # the host order is inherited: keep is already in cyclic order
    pos = {x: i for i, x in enumerate(keep)}
    m2 = [-1]*len(keep)
    for x in keep:
        y = mate[x]
        if y in dead: return None      # cannot happen: mate[x] != u,v for x kept
        m2[pos[x]] = pos[y]
    return m2

def girth(mate, n):
    A = [[(v+1) % n, (v-1) % n, mate[v]] for v in range(n)]
    g = 10**9
    for src in range(n):
        d = {src: 0}; par = {src: -1}; q = deque([src])
        while q:
            x = q.popleft()
            if 2*d[x] >= g: break
            for w in A[x]:
                if w not in d:
                    d[w] = d[x]+1; par[w] = x; q.append(w)
                elif w != par[x]: g = min(g, d[x]+d[w]+1)
    return g

def valid_ladder(mate, n):
    if n < 6: return False
    if any(mate[mate[x]] != x or mate[x] == x for x in range(n)): return False
    if any(mate[x] in ((x+1) % n, (x-1) % n) for x in range(n)): return False
    return True

if __name__ == '__main__':
    from skip3 import random_ladder, all_ladders
    rng = random.Random(1618)
    print("contraction: how many chords give a girth-legal smaller ladder?")
    print("  2N   ladders   min good chords   ladders with NONE")
    for n in (12, 14, 16, 18):
        Ls = [m for m in all_ladders(n) if girth(m, n) >= 5]
        worst = 10**9; none = 0
        for m in Ls:
            good = 0
            for u in range(n):
                if mate_is_left := (u < m[u]):
                    m2 = contract(m, n, u)
                    if m2 is not None and valid_ladder(m2, n-2) \
                       and girth(m2, n-2) >= 5:
                        good += 1
            worst = min(worst, good)
            if good == 0: none += 1
        print("  %2d   %7d        %3d             %d" % (n, len(Ls), worst, none))
        sys.stdout.flush()
