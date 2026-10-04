# -*- coding: utf-8 -*-
"""516: verify the 2N = 20 residue with the independent template
recogniser of 515 (k33b.contract / k33b.is_target), which knows
nothing about crossings, spans or weight formulas.

For each closed girth-5 ladder with k* = 3, and for each three-chord
power-of-two cycle Z in it, ask: is C + (the three chords Z uses) a
K_{3,3}-subdivision?  If not, then Z lies in NO template through C -
no theta, no K_4, no K_{3,3} - and every through-C statement in the
archive is silent on the only power-of-two cycle the ladder has.
"""

import sys
from itertools import combinations
from kstar2 import cyc_lengths, sweep
from k33b import contract, is_target


def host_plus(chords, n):
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    for a, b in chords:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def is_k33_sub(chords, n):
    adj = host_plus(chords, n)
    br, ed = contract(adj, sorted(adj))
    return is_target(br, ed, 6)


def is_k4_sub(chords, n):
    adj = host_plus(chords, n)
    br, ed = contract(adj, sorted(adj))
    return is_target(br, ed, 4)


def audit(n):
    surv, hist, wit = sweep(n, show=0)
    P2 = set()
    v = 4
    while v <= 2 * n:
        P2.add(v)
        v *= 2
    covered = uncovered = 0
    both = 0
    for ch in surv:
        temp_hit = False
        any_hit = False
        for T in combinations(ch, 3):
            got = cyc_lengths(list(T), n, None)
            if any(c == 3 and L in P2 for L, c in got.items()):
                any_hit = True
                if is_k33_sub(list(T), n):
                    temp_hit = True
        assert any_hit
        if temp_hit:
            covered += 1
        else:
            uncovered += 1
    # sanity: on ALL triples of ALL survivors, crossing <=> K_{3,3}
    mism = 0
    for ch in surv[:20]:
        for T in combinations(ch, 3):
            (a1, b1), (a2, b2), (a3, b3) = sorted(T)
            cross = a1 < a2 < a3 < b1 < b2 < b3
            if cross != is_k33_sub(list(T), n):
                mism += 1
    print("2N=%2d  k*=3 ladders=%4d   power of two inside a through-C "
          "K_{3,3}: %4d   in NO through-C template: %4d"
          % (n, len(surv), covered, uncovered))
    print("       cross-check on 20 ladders, all triples: "
          "'pairwise crossing' vs 'is a K_{3,3}-subdivision' "
          "mismatches = %d" % mism)
    sys.stdout.flush()
    return uncovered


if __name__ == '__main__':
    audit(18)
    audit(20)
