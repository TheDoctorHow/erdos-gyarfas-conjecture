# -*- coding: utf-8 -*-
"""Machine audit for section 500: the zigzag calculus on girth-legal
ladders.  Everything is recomputed from the edge set with independent
code (union-find, explicit degree check); no formula is trusted."""

import random
from itertools import combinations

# ---------------------------------------------------------------- model

def overlap(c, d):
    (a, b), (x, y) = c, d
    return not (b < x or y < a)

def legal_pair(c, d):
    if not overlap(c, d):
        return True
    (a, b), (x, y) = c, d
    return abs(a - x) + abs(b - y) >= 3

def gen(N, rng, tries=4000):
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

# ------------------------------------------------- completion, honestly

def completion(T, n):
    """Edges of E(T): the chords of T plus the unique host edges making
    every degree even.  Returns (edges, weight, connected)."""
    ends = sorted([p for c in T for p in c])
    edges = [tuple(c) for c in T]
    cnt = 0
    idx = 0
    for j in range(1, n):                     # host edge (j, j+1)
        while idx < len(ends) and ends[idx] <= j:
            idx += 1
            cnt += 1
        if cnt % 2 == 1:
            edges.append((j, j + 1))
    deg = {}
    for u, v in edges:
        deg[u] = deg.get(u, 0) + 1
        deg[v] = deg.get(v, 0) + 1
    assert all(d % 2 == 0 for d in deg.values()), "completion not even"
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
    roots = set(find(v) for v in deg)
    conn = (len(roots) <= 1) and len(deg) > 0
    return edges, len(edges), conn

def alt_weight(T):
    q = sorted([p for c in T for p in c])
    return len(T) + sum((-1) ** (i + 1) * q[i] for i in range(len(q)))

# --------------------------------------------------------------- claims

def eps(c):
    a, b = c
    return 1 + (a % 2) - (b % 2)

def straddle(chords, j):
    return sum(1 for (a, b) in chords if a <= j < b)

def greedy_chain(chords, n):
    """Greedy zigzag chain: start at the chord of position 1, then
    repeatedly take the straddling chord of largest right end."""
    cur = [c for c in chords if c[0] == 1][0]
    chain = [cur]
    while chain[-1][1] < n:
        b = chain[-1][1]
        cands = [c for c in chords if c[0] <= b < c[1]]
        if not cands:
            return None                       # reducible ladder
        chain.append(max(cands, key=lambda c: c[1]))
    return chain

def sub_weight(chain, i, j):
    """W(i,j) by the (count) + (b_j - a_i) - sum overlaps formula."""
    cnt = j - i + 1
    span = chain[j][1] - chain[i][0]
    ov = sum(chain[t][1] - chain[t + 1][0] for t in range(i, j))
    return cnt + span - ov

def RL(chain):
    O = [0]
    for t in range(len(chain) - 1):
        O.append(O[-1] + chain[t][1] - chain[t + 1][0])
    R = [(j + 1) + 1 + chain[j][1] - O[j] for j in range(len(chain))]
    L = [(i + 1) + chain[i][0] - O[i] for i in range(len(chain))]
    return R, L

# ----------------------------------------------------------------- test

def main():
    rng = random.Random(20260809)
    stats = dict((k, 0) for k in
                 ['ladders', 'wform', 'eps', 'disj', 'parity', 'toggle',
                  'strad', 'chains', 'zz', 'sub', 'rl', 'reach',
                  'swap', 'ins'])
    fails = []

    for N in range(5, 12):
        for _ in range(40):
            ch = gen(N, rng)
            if ch is None:
                continue
            n = 2 * N
            stats['ladders'] += 1

            # 1. weight formula, on random subsets
            for _ in range(30):
                k = rng.randint(1, N)
                T = rng.sample(ch, k)
                _, w, _ = completion(T, n)
                if w != alt_weight(T):
                    fails.append(('wform', ch, T))
                else:
                    stats['wform'] += 1
                # parity invariant
                p = (len(T) + sum(b - a for a, b in T)) % 2
                if w % 2 != p:
                    fails.append(('parity', ch, T))
                else:
                    stats['parity'] += 1

            # 2. epsilon formula for single-chord deletion
            for c in ch:
                T = [d for d in ch if d != c]
                _, w, _ = completion(T, n)
                if w != n - eps(c):
                    fails.append(('eps', ch, c))
                else:
                    stats['eps'] += 1

            # 3. disjoint-family additivity + reach limit
            best = 0
            for r in range(1, 5):
                for S in combinations(ch, r):
                    if any(overlap(x, y) for x, y in combinations(S, 2)):
                        continue
                    T = [d for d in ch if d not in S]
                    _, w, _ = completion(T, n)
                    if w != n - sum(eps(c) for c in S):
                        fails.append(('disj', ch, S))
                    else:
                        stats['disj'] += 1
                    best = max(best, sum(eps(c) for c in S))
            if best > 2 * (n - 1) / 5.0:
                fails.append(('reach', ch, best))
            else:
                stats['reach'] += 1

            # 4. toggle formula  dw = 1 + sigma - 2u
            for _ in range(30):
                k = rng.randint(0, N - 1)
                T = rng.sample(ch, k)
                rest = [d for d in ch if d not in T]
                c = rng.choice(rest)
                e0, w0, _ = completion(T, n)
                _, w1, _ = completion(T + [c], n)
                used = set(x for x in e0 if x[1] == x[0] + 1)
                u = sum(1 for j in range(c[0], c[1]) if (j, j + 1) in used)
                if w1 - w0 != 1 + (c[1] - c[0]) - 2 * u:
                    fails.append(('toggle', ch, T, c))
                else:
                    stats['toggle'] += 1

            # 5. straddle parity
            if all(straddle(ch, j) % 2 == j % 2 for j in range(1, n)):
                stats['strad'] += 1
            else:
                fails.append(('strad', ch))

            # 6. the zigzag chain
            chain = greedy_chain(ch, n)
            if chain is None:
                continue
            stats['chains'] += 1
            m = len(chain)
            good = True
            for t in range(m - 1):
                a, b = chain[t]
                x, y = chain[t + 1]
                if not (a < x < b < y):
                    good = False
            for t in range(m - 2):
                if not (chain[t][1] < chain[t + 2][0]):
                    good = False
            if good:
                stats['zz'] += 1
            else:
                fails.append(('zz', ch, chain))

            R, L = RL(chain)
            for i in range(m):
                for j in range(i, m):
                    T = chain[i:j + 1]
                    _, w, conn = completion(T, n)
                    if not conn:
                        fails.append(('sub', ch, i, j))
                    else:
                        stats['sub'] += 1
                    if w != sub_weight(chain, i, j) or w != R[j] - L[i]:
                        fails.append(('rl', ch, i, j, w))
                    else:
                        stats['rl'] += 1

            # 7. the three moves: insertion and interior swap
            for t in range(m - 1):
                a, b = chain[t]
                x, y = chain[t + 1]
                for c in ch:
                    if c in chain:
                        continue
                    p, q = c
                    if not (a < p < b < q < y and p > (chain[t - 1][1]
                            if t > 0 else 0)
                            and (t + 2 >= m or q < chain[t + 2][0])):
                        continue
                    new = chain[:t + 1] + [c] + chain[t + 1:]
                    _, w0, _ = completion(chain, n)
                    _, w1, cc = completion(new, n)
                    if cc and w1 - w0 != 1 - (q - p):
                        fails.append(('ins', ch, t, c))
                    elif cc:
                        stats['ins'] += 1
            for t in range(1, m - 1):
                a0, b0 = chain[t - 1]
                a2, b2 = chain[t + 1]
                for c in ch:
                    if c in chain:
                        continue
                    p, q = c
                    if not (a0 < p < b0 < q and p < a2 < q < b2):
                        continue
                    if t - 2 >= 0 and not chain[t - 2][1] < p:
                        continue
                    if t + 2 < m and not q < chain[t + 2][0]:
                        continue
                    new = chain[:t] + [c] + chain[t + 1:]
                    _, w0, _ = completion(chain, n)
                    _, w1, cc = completion(new, n)
                    d = (chain[t][1] - chain[t][0]) - (q - p)
                    if cc and w1 - w0 != d:
                        fails.append(('swap', ch, t, c))
                    elif cc:
                        stats['swap'] += 1

    print("stats:", stats)
    print("failures:", len(fails))
    for f in fails[:8]:
        print("  ", f)

main()
