# -*- coding: utf-8 -*-
"""Exhaustive search for a (K4) counterexample.

"Some K4-subdivision has a power of two among its seven" is MONOTONE
in the chord set - adding chords never destroys a subdivision - so its
negation prunes soundly.  We prune with a CHEAP certificate (the
crossing-triple K4s, whose seven are the one-chord, two-chord and XXX
lengths of 508.2) and then run the FULL test on every survivor.
"""

import sys
from itertools import combinations

POW = set(1 << i for i in range(1, 24))
exec(open('k4c.py').read().split('def main()')[0].split(
    '# -*- coding: utf-8 -*-')[1].split('"""', 2)[2])

def cheap_hit(chords, newc):
    """Is there a crossing triple / pair / single giving a power of
    two, involving newc?  (== the seven of a short-branch K4.)"""
    a, b = newc
    if (b - a + 1) in POW:
        return True
    for d in chords:
        if d == newc:
            continue
        if ov(newc, d) and (abs(a - d[0]) + abs(b - d[1]) + 2) in POW:
            return True
    others = [c for c in chords if c != newc]
    for c1, c2 in combinations(others, 2):
        for T in ((c1, c2, newc), (c1, newc, c2), (newc, c1, c2)):
            (a1, b1), (a2, b2), (a3, b3) = T
            if a1 < a2 < a3 < b1 < b2 < b3:
                if (3 + (b1 - a1) - (b2 - a2) + (b3 - a3)) in POW:
                    return True
    return False

def search(n):
    used = [False] * (n + 2)
    chords = []
    nodes = [0]
    survivors = [0]
    found = []

    def rec():
        nodes[0] += 1
        p = 1
        while p <= n and used[p]:
            p += 1
        if p > n:
            survivors[0] += 1
            ok, _ = k4_holds(list(chords), n)
            if not ok:
                found.append(list(chords))
                return True
            return False
        used[p] = True
        for q in range(p + 4, n + 1):
            if used[q]:
                continue
            c = (p, q)
            if any(not (abs(c[0] - d[0]) + abs(c[1] - d[1]) >= 3)
                   for d in chords if ov(c, d)):
                continue
            used[q] = True
            chords.append(c)
            if not cheap_hit(chords, c):
                if rec():
                    return True
            chords.pop()
            used[q] = False
        used[p] = False
        return False

    hit = rec()
    return (found[0] if hit else None), nodes[0], survivors[0]

def main():
    sys.setrecursionlimit(30000)
    for n in range(12, 31, 2):
        w, nodes, surv = search(n)
        if w is None:
            print("2N =", n, " no (K4) counterexample; nodes", nodes,
                  " survivors fully tested:", surv)
        else:
            print("2N =", n, " *** (K4) COUNTEREXAMPLE ***", w,
                  " nodes", nodes)
            break

main()
