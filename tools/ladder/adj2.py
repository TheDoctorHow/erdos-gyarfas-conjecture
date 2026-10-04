# -*- coding: utf-8 -*-
"""Adversarial search: minimise the number of (ADJ) witnesses over
hard-case ladders.  If a hill-climb can drive the witness count to 0,
(ADJ) is false; if the floor stays well above 0 and grows with 2N,
(ADJ) is robust and we learn what the extremal ladders look like.

Move set: pick two chords and re-pair their four endpoints (the two
alternative pairings), keep if still hard-case and the witness count
does not increase.
"""

import sys, random

POW = set(1 << i for i in range(1, 22))

def overlap(c, d):
    (a, b), (p, q) = c, d
    return not (b < p or q < a)

def l1(c, d):
    return abs(c[0] - d[0]) + abs(c[1] - d[1])

def hard_pair(c, d):
    if not overlap(c, d):
        return True
    v = l1(c, d)
    return v >= 3 and (v + 2) not in POW

def hard_span(c):
    s = c[1] - c[0]
    return s >= 4 and (s + 1) not in POW

def legal(ch):
    for i in range(len(ch)):
        if not hard_span(ch[i]):
            return False
        for j in range(i + 1, len(ch)):
            if not hard_pair(ch[i], ch[j]):
                return False
    return True

def witnesses(ch, n):
    mu = {}
    for a, b in ch:
        mu[a] = b
        mu[b] = a
    cnt = 0
    for p in range(1, n):
        c = tuple(sorted((p, mu[p])))
        cp = tuple(sorted((p + 1, mu[p + 1])))
        if c == cp:
            continue
        Q = (min(mu[p], mu[p + 1]), max(mu[p], mu[p + 1]))
        beta = (mu[p] < p) + (mu[p + 1] < p)
        for D in ch:
            if D == c or D == cp:
                continue
            if not overlap(D, Q):
                continue
            if (beta + (D[0] < p) + (D[1] < p)) % 2:
                continue
            if 4 + l1(Q, D) in POW:
                cnt += 1
    return cnt

def gen_hard(n, rng, tries=20000):
    for _ in range(tries):
        used = [False] * (n + 2)
        chords = []
        ok = True
        for p in range(1, n + 1):
            if used[p]:
                continue
            used[p] = True
            cand = [q for q in range(p + 4, n + 1) if not used[q]]
            rng.shuffle(cand)
            placed = False
            for q in cand:
                c = (p, q)
                if not hard_span(c):
                    continue
                if any(not hard_pair(c, d) for d in chords):
                    continue
                used[q] = True
                chords.append(c)
                placed = True
                break
            if not placed:
                ok = False
                break
        if ok and len(chords) * 2 == n:
            return sorted(chords)
    return None

def repair(ch, i, j, mode):
    (a, b), (c, d) = ch[i], ch[j]
    pts = sorted([a, b, c, d])
    if mode == 0:
        new = [(pts[0], pts[1]), (pts[2], pts[3])]
    elif mode == 1:
        new = [(pts[0], pts[2]), (pts[1], pts[3])]
    else:
        new = [(pts[0], pts[3]), (pts[1], pts[2])]
    out = [ch[k] for k in range(len(ch)) if k not in (i, j)] + new
    return sorted(out)

def climb(ch, n, rng, steps=40000):
    best = witnesses(ch, n)
    cur = ch
    for _ in range(steps):
        i, j = rng.sample(range(len(cur)), 2)
        m = rng.randint(0, 2)
        nxt = repair(cur, i, j, m)
        if len(set(x for c in nxt for x in c)) != n:
            continue
        if not legal(nxt):
            continue
        v = witnesses(nxt, n)
        if v <= best:
            if v < best or rng.random() < 0.35:
                cur, best = nxt, v
        if best == 0:
            break
    return cur, best

def main():
    rng = random.Random(24680)
    for n in (20, 26, 32, 40, 50):
        floor = None
        arg = None
        for r in range(6):
            ch = gen_hard(n, rng)
            if ch is None:
                continue
            c2, v = climb(ch, n, rng)
            if floor is None or v < floor:
                floor, arg = v, c2
            if v == 0:
                break
        print("2N =", n, " minimum (ADJ) witness count found:", floor)
        print("    extremal ladder:", arg)
        if floor == 0:
            print("    *** (ADJ) IS FALSE ***")

main()
