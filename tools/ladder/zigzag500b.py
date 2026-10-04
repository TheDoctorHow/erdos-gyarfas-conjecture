# -*- coding: utf-8 -*-
"""Second audit for section 500.

(1) STRADDLE-PARITY FORMULA: host edge (j,j+1) is in the completion of
    T iff an odd number of T-chords straddle gap j; hence
    w(T) = |T| + #{j : k_T(j) odd}.
(2) STAIRCASE HYPOTHESIS: if the chords of T can be listed with both
    endpoint sequences strictly increasing and consecutive ones
    crossing, is E(T) connected?
(3) Do staircases alone realise a power-of-two cycle in every
    girth-legal ladder?
"""

import random
from itertools import combinations

def overlap(c, d):
    (a, b), (x, y) = c, d
    return not (b < x or y < a)

def legal_pair(c, d):
    if not overlap(c, d):
        return True
    (a, b), (x, y) = c, d
    return abs(a - x) + abs(b - y) >= 3

def gen(N, rng, tries=6000):
    n = 2 * N
    for _ in range(tries):
        partner, chords, ok = {}, [], True
        for p in range(1, n + 1):
            if p in partner:
                continue
            cands = [q for q in range(p + 4, n + 1) if q not in partner]
            rng.shuffle(cands)
            placed = False
            for q in cands:
                c = (p, q)
                if all(legal_pair(c, d) for d in chords):
                    partner[p] = q
                    partner[q] = p
                    chords.append(c)
                    placed = True
                    break
            if not placed:
                ok = False
                break
        if ok and len(chords) == N:
            return sorted(chords)
    return None

def completion(T, n):
    ends = sorted([p for c in T for p in c])
    edges = [tuple(c) for c in T]
    cnt, idx = 0, 0
    for j in range(1, n):
        while idx < len(ends) and ends[idx] <= j:
            idx += 1
            cnt += 1
        if cnt % 2 == 1:
            edges.append((j, j + 1))
    deg = {}
    for u, v in edges:
        deg[u] = deg.get(u, 0) + 1
        deg[v] = deg.get(v, 0) + 1
    assert all(d % 2 == 0 for d in deg.values())
    par = {v: v for v in deg}

    def find(v):
        while par[v] != v:
            par[v] = par[par[v]]
            v = par[v]
        return v

    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            par[ru] = rv
    conn = len(set(find(v) for v in deg)) <= 1 and len(deg) > 0
    return len(edges), conn

def straddle_weight(T, n):
    return len(T) + sum(1 for j in range(1, n)
                        if sum(1 for (a, b) in T if a <= j < b) % 2 == 1)

def is_staircase(T):
    S = sorted(T)
    for t in range(len(S) - 1):
        if not (S[t][0] < S[t + 1][0] and S[t][1] < S[t + 1][1]):
            return False
        if not S[t + 1][0] < S[t][1]:
            return False
    return True

def pow2(x):
    return x >= 2 and (x & (x - 1)) == 0

def main():
    rng = random.Random(4711)
    nform = nstair = nstair_ok = 0
    stair_fail = []
    lad = 0
    all_p2 = stair_p2 = 0
    no_p2_at_all = []
    covered = 0
    for N in range(5, 11):
        for _ in range(60):
            ch = gen(N, rng)
            if ch is None:
                continue
            n = 2 * N
            lad += 1
            specA, specS = set(), set()
            for r in range(1, N + 1):
                for T in combinations(ch, r):
                    w, conn = completion(list(T), n)
                    if w != straddle_weight(list(T), n):
                        print("STRADDLE FORMULA FAILS", ch, T)
                    else:
                        nform += 1
                    st = is_staircase(T)
                    if st:
                        nstair += 1
                        if conn:
                            nstair_ok += 1
                        else:
                            stair_fail.append((ch, T))
                    if conn:
                        specA.add(w)
                        if st:
                            specS.add(w)
            a = any(pow2(x) for x in specA)
            s = any(pow2(x) for x in specS)
            all_p2 += a
            stair_p2 += s
            if not a:
                no_p2_at_all.append(ch)
            if specS >= specA:
                covered += 1
    print("ladders", lad)
    print("straddle-formula checks", nform)
    print("staircases", nstair, "connected", nstair_ok,
          "fails", len(stair_fail))
    for f in stair_fail[:5]:
        print("   ", f)
    print("ladders with a power-of-two cycle:", all_p2)
    print("  ... realised by a staircase:", stair_p2)
    print("ladders whose staircase spectrum = full spectrum:", covered)
    for f in no_p2_at_all[:5]:
        print("  NO POWER OF TWO:", f)

main()
