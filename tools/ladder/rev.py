import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
# Reverse-engineering check.
#
# CLAIM (the fully reduced target).  Let H be the cubic graph of the normal form
# 224.2 and e1,e2,e3 the matching whose subdivision points are the terminals.
# Then for the pair (b_i,b_j):
#
#    two internally disjoint EVEN connecting paths avoiding b_k
#      <==>  H has a cycle through e_i and e_j, NOT through e_k,
#            whose two e_i-e_j arcs are BOTH EVEN.
#
# Verified here both ways: the H-side by cycle enumeration, the B-side by brute
# force over connecting paths in the subdivided graph.
import random, itertools, sys
from collections import deque

base = _HERE
src = open(base + "surv.py").read()
exec(src[:src.rindex("def main")])


def all_cycles(adj, n):
    out = []
    for s in range(n):
        stack = [(s, [s], {s})]
        while stack:
            u, path, seen = stack.pop()
            for w in adj[u]:
                if w == s and len(path) >= 3:
                    if path[1] < path[-1]:
                        out.append(tuple(path))
                elif w > s and w not in seen:
                    stack.append((w, path + [w], seen | {w}))
    return out


def cyc_edges(cyc):
    k = len(cyc)
    return [frozenset((cyc[i], cyc[(i + 1) % k])) for i in range(k)]


def h_side(cycles, e1, e2, e3):
    """cycle through e1,e2, avoiding e3, both e1-e2 arcs even?"""
    f1, f2, f3 = frozenset(e1), frozenset(e2), frozenset(e3)
    for cyc in cycles:
        E = cyc_edges(cyc)
        Es = set(E)
        if f1 not in Es or f2 not in Es or f3 in Es:
            continue
        k = len(cyc)
        i1 = E.index(f1)
        i2 = E.index(f2)
        a = (i2 - i1 - 1) % k          # edges strictly between e1 and e2
        b = k - 2 - a
        if a % 2 == 0 and b % 2 == 0:
            return True
    return False


def build_B(nH, adj, M):
    """subdivide the three matching edges; terminals are nH, nH+1, nH+2"""
    badj = [list(x) for x in adj]
    for i, (a, b) in enumerate(M):
        t = nH + i
        badj[a].remove(b)
        badj[b].remove(a)
        badj.append([a, b])
        badj[a].append(t)
        badj[b].append(t)
    return badj


def b_side(badj, N, b1, b2, b3):
    """Independent check on the SUBDIVIDED graph B: is there a cycle through
    b1 and b2, avoiding b3, with both b1-b2 arcs EVEN?  Enumerated by walking
    from one neighbour of b1 round to the other."""
    p, q = badj[b1]
    stack = [(p, (1 << b1) | (1 << p), [p])]
    while stack:
        u, vis, path = stack.pop()
        for w in badj[u]:
            if w == q and len(path) >= 2:
                cyc = [b1] + path + [q]
                if b2 in cyc:
                    k = len(cyc)
                    i2 = cyc.index(b2)
                    if i2 % 2 == 0 and (k - i2) % 2 == 0:
                        return True
                continue
            if w == b3 or w == b1 or (vis >> w) & 1:
                continue
            stack.append((w, vis | (1 << w), path + [w]))
    return False


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
    mismatch = 0
    cfg = 0
    from collections import Counter
    per = Counter()
    zerowork = 0
    checked_B = 0
    for nH in (14, 16, 18, 20):
        for trial in range(60):
            g = randcubic(nH, rng)
            if not g:
                continue
            E, adj = g
            D = [dist(adj, nH, v) for v in range(nH)]

            def edist(e, f):
                return min(D[a][b] for a in e for b in f)
            cand = [t for t in itertools.combinations(E, 3)
                    if all(edist(t[i], t[j]) >= 3
                           for i in range(3) for j in range(i + 1, 3))]
            if not cand:
                continue
            rng.shuffle(cand)
            cycles = all_cycles(adj, nH)
            for M in cand[:3]:
                cfg += 1
                works = 0
                for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
                    hs = h_side(cycles, M[i], M[j], M[k])
                    if hs:
                        works += 1
                    badj = build_B(nH, adj, M)
                    bs = b_side(badj, nH + 3, nH + i, nH + j, nH + k)
                    checked_B += 1
                    if bs != hs:
                        mismatch += 1
                per[works] += 1
                if works == 0:
                    zerowork += 1
    print("configurations (H, matching of 3, pairwise dist >= 3) :", cfg)
    print("pairs cross-checked H-side vs B-side                  :", checked_B)
    print("  mismatches                                          :", mismatch)
    print()
    print("number of the three pairs admitting the target cycle:")
    for k in sorted(per):
        print("   %d pairs work -> %5d configurations" % (k, per[k]))
    print("configurations with NO working pair                   :", zerowork)


main()
