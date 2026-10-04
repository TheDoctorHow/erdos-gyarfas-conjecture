"""(WIND0-GEN) at the top of the range.

Taking ALL rungs of one period gives coverage parity chi(t) = t+1 mod 2,
so every other vertex is used and the length is 2N.  REMOVING a rung
INVERTS the parity on its span, so with T the removed set

    length  =  2N  -  |T|  +  D(T),
    D(T)    =  #{t odd, covered oddly by T} - #{t even, covered oddly}.

The odd-covered set of T is exactly T's own (WIND0) pairing intervals,
so D is a sum of one term per interval and the DEFICIENCY is

    2N - length  =  sum over T's pairing intervals of c_i,
    c_i = 1 if the interval has even length,
          0 if odd length starting at an odd position,
          2 if odd length starting at an even position.
"""
import sys, random, itertools
from gwind import lift, girth

def rungs(mate, n):
    isleft, span = lift(mate, n)
    return [(x, span[x]) for x in range(n) if isleft[x]]

def wind0_length(S, n):
    """S = list of (x, sigma) in the cover.  Returns (length, ok, verts)."""
    k = len(S)
    pts = []
    for i, (x, s) in enumerate(S):
        pts.append((x, i)); pts.append((x+s, i))
    pts.sort()
    if len(set(p for p, _ in pts)) != 2*k: return None, False, None
    adj = [[] for _ in range(k)]
    host = 0; verts = set()
    for j in range(0, 2*k, 2):
        (a, ia), (b, ib) = pts[j], pts[j+1]
        adj[ia].append(ib); adj[ib].append(ia)
        host += b-a
        for t in range(a, b+1): verts.add(t)
    if any(len(x) != 2 for x in adj): return None, False, None
    seen = {0}; prev = -1; cur = 0
    for _ in range(k):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        prev, cur = cur, nxt; seen.add(cur)
    return k+host, len(seen) == k, verts

def deficiency_formula(T, n):
    """sum of c_i over T's pairing intervals."""
    if not T: return 0
    pts = sorted([x for x, s in T] + [x+s for x, s in T])
    tot = 0
    for j in range(0, len(pts), 2):
        a, b = pts[j], pts[j+1]
        L = b-a
        if L % 2 == 0: tot += 1
        elif a % 2 == 1: tot += 0
        else: tot += 2
    return tot

if __name__ == '__main__':
    from skip3 import random_ladder
    rng = random.Random(60221)
    print("== the dual formula: length = 2N - deficiency(T) ==")
    bad = 0; tested = 0
    for n in (16, 20, 24, 28):
        for _ in range(4):
            m = random_ladder(n, rng)
            if m is None or girth(m, n) < 5: continue
            R = rungs(m, n)
            for tsize in range(0, min(4, len(R))+1):
                for T in itertools.combinations(range(len(R)), tsize):
                    S = [R[i] for i in range(len(R)) if i not in T]
                    Tl = [R[i] for i in T]
                    L, ok, verts = wind0_length(S, n)
                    if L is None: continue
                    pred = n - deficiency_formula(Tl, n)
                    tested += 1
                    if L != pred:
                        bad += 1
                        if bad <= 4:
                            print("   MISMATCH n=%d |T|=%d  actual %d predicted %d"
                                  % (n, tsize, L, pred))
    print("  configurations tested %d   formula mismatches %d" % (tested, bad))
