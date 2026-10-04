# -*- coding: utf-8 -*-
"""Third audit: is the ZIGZAG family alone enough for (LQ3)?

A zigzag family is a chord set that, sorted by left endpoint, has both
endpoint sequences increasing, consecutive pairs crossing, and
non-consecutive pairs DISJOINT.  Claim A: every zigzag family has a
connected completion.  Question B: does every girth-legal ladder have a
zigzag family of power-of-two weight?
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

def is_zigzag(T):
    S = sorted(T)
    m = len(S)
    for t in range(m - 1):
        if not (S[t][0] < S[t + 1][0] < S[t][1] < S[t + 1][1]):
            return False
    for t in range(m - 2):
        if not S[t][1] < S[t + 2][0]:
            return False
    return True

def pow2(x):
    return x >= 2 and (x & (x - 1)) == 0

def main():
    rng = random.Random(99001)
    lad = 0
    zz = zz_conn = 0
    zz_bad = []
    have_all = have_zz = 0
    misses = []
    for N in range(5, 12):
        for _ in range(60):
            ch = gen(N, rng)
            if ch is None:
                continue
            n = 2 * N
            lad += 1
            specA, specZ = set(), set()
            for r in range(1, N + 1):
                for T in combinations(ch, r):
                    w, conn = completion(list(T), n)
                    if conn:
                        specA.add(w)
                    if is_zigzag(T):
                        zz += 1
                        if conn:
                            zz_conn += 1
                            specZ.add(w)
                        else:
                            zz_bad.append((ch, T))
            a = any(pow2(x) for x in specA)
            z = any(pow2(x) for x in specZ)
            have_all += a
            have_zz += z
            if a and not z:
                misses.append((ch, sorted(specA), sorted(specZ)))
    print("ladders", lad)
    print("zigzag families", zz, "connected", zz_conn,
          "disconnected", len(zz_bad))
    for b in zz_bad[:5]:
        print("   BAD", b)
    print("ladders with a power-of-two cycle:", have_all)
    print("  ... realised by a zigzag family:", have_zz)
    print("  misses:", len(misses))
    for m in misses[:6]:
        print("   MISS", m[0])
        print("        all:", m[1])
        print("        zig:", m[2])

main()
