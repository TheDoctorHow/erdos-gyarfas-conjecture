# -*- coding: utf-8 -*-
"""Can the planar theorem reach a ladder?

A ladder = host path + perfect matching has a 2-page book embedding
(hence is planar) as soon as the CROSSING GRAPH of its chords is
bipartite: colour the chords by page.  So:

 (1) how many girth-legal ladders are planar in this sense?
 (2) the dodecahedron minus a Hamiltonian-cycle edge should BE one -
     build it from LCF notation, verify it is girth-legal, and read
     off its spectrum.
"""

import sys
from itertools import combinations
from collections import deque

def overlap(c, d):
    (a, b), (p, q) = c, d
    return not (b < p or q < a)

def crosses(c, d):
    (a, b), (p, q) = c, d
    return a < p < b < q or p < a < q < b

def soft_pair(c, d):
    if not overlap(c, d):
        return True
    return abs(c[0] - d[0]) + abs(c[1] - d[1]) >= 3

def bipartite_crossing(ch):
    n = len(ch)
    col = [-1] * n
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if crosses(ch[i], ch[j]):
                adj[i].append(j)
                adj[j].append(i)
    for s in range(n):
        if col[s] != -1:
            continue
        col[s] = 0
        Q = deque([s])
        while Q:
            v = Q.popleft()
            for w in adj[v]:
                if col[w] == -1:
                    col[w] = 1 - col[v]
                    Q.append(w)
                elif col[w] == col[v]:
                    return False
    return True

def enumerate_soft(n):
    out, chords = [], []
    used = [False] * (n + 2)

    def rec():
        p = 1
        while p <= n and used[p]:
            p += 1
        if p > n:
            out.append(list(chords))
            return
        used[p] = True
        for q in range(p + 4, n + 1):
            if used[q]:
                continue
            c = (p, q)
            if any(not soft_pair(c, d) for d in chords):
                continue
            used[q] = True
            chords.append(c)
            rec()
            chords.pop()
            used[q] = False
        used[p] = False

    rec()
    return out

def completion(T, n):
    ends = sorted([x for c in T for x in c])
    edges = [tuple(c) for c in T]
    cnt, idx = 0, 0
    for j in range(1, n):
        while idx < len(ends) and ends[idx] <= j:
            idx += 1
            cnt += 1
        if cnt % 2 == 1:
            edges.append((j, j + 1))
    deg = {}
    for a, b in edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    par = {v: v for v in deg}

    def find(v):
        while par[v] != v:
            par[v] = par[par[v]]
            v = par[v]
        return v

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb
    return len(edges), len(set(find(v) for v in deg)) <= 1

def spectrum(ch, n):
    S = set()
    for r in range(1, len(ch) + 1):
        for T in combinations(ch, r):
            w, c = completion(list(T), n)
            if c:
                S.add(w)
    return sorted(S)

# ---------------------------------------------------- the dodecahedron

LCF = [10, 7, 4, -4, -7, 10, -4, 7, -7, 4] * 2

def dodeca():
    n = 20
    E = set()
    for i in range(n):
        E.add(tuple(sorted((i, (i + 1) % n))))
    for i, d in enumerate(LCF):
        E.add(tuple(sorted((i, (i + d) % n))))
    return n, sorted(E)

def girth_of(n, E):
    adj = {v: [] for v in range(n)}
    for a, b in E:
        adj[a].append(b)
        adj[b].append(a)
    best = 10 ** 9
    for s in range(n):
        dist = {s: 0}
        par = {s: None}
        Q = deque([s])
        while Q:
            v = Q.popleft()
            for w in adj[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1
                    par[w] = v
                    Q.append(w)
                elif par[v] != w:
                    best = min(best, dist[v] + dist[w] + 1)
    return best

def main():
    sys.setrecursionlimit(20000)
    print("=== (1) planar girth-legal ladders")
    for n in (12, 14, 16):
        lads = enumerate_soft(n)
        pl = [ch for ch in lads if bipartite_crossing(ch)]
        print("2N =", n, " girth-legal", len(lads),
              " bipartite crossing graph (=> planar):", len(pl))
        for ch in pl[:2]:
            print("      ", ch, " spectrum", spectrum(ch, n))

    print()
    print("=== (2) the dodecahedron minus a Hamiltonian-cycle edge")
    n, E = dodeca()
    deg = {v: 0 for v in range(n)}
    for a, b in E:
        deg[a] += 1
        deg[b] += 1
    print("vertices", n, " edges", len(E),
          " all degree 3?", set(deg.values()) == {3},
          " girth", girth_of(n, E))
    # remove the cycle edge (19,0); relabel i -> i+1
    chords = []
    for a, b in E:
        if (a, b) == (0, 19):
            continue
        if b == a + 1:
            continue
        chords.append(tuple(sorted((a + 1, b + 1))))
    chords = sorted(chords)
    print("as a ladder on [1,20], chords:", chords)
    print("spans:", sorted(b - a for a, b in chords))
    viol = [(c, d) for c, d in combinations(chords, 2)
            if not soft_pair(c, d)]
    print("girth-legal?", all(b - a >= 4 for a, b in chords)
          and not viol, " violations", viol[:3])
    print("crossing graph bipartite (=> planar)?",
          bipartite_crossing(chords))
    S = spectrum(chords, 20)
    print("spectrum:", S)
    print("girth", S[0], " circumference", S[-1],
          "  8 present?", 8 in S, "  16 present?", 16 in S)

main()
