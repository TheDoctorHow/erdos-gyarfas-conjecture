import itertools, sys
from math import gcd
from collections import deque

# ---------- matchings ----------
def matchings(verts):
    if not verts:
        yield []
        return
    a = verts[0]
    for i in range(1, len(verts)):
        b = verts[i]
        rest = verts[1:i] + verts[i+1:]
        for m in matchings(rest):
            yield [(a, b)] + m

def rel(e, f):
    a, b = min(e), max(e)
    c, d = min(f), max(f)
    if a < c < b < d or c < a < d < b:
        return 'X'
    if a < c < d < b or c < a < b < d:
        return 'N'
    return 'D'

def open_legal(m):
    """spans >= 4 and non-disjoint pairs at L1 >= 3: exactly the
    conditions a self-matched arc of a girth-5 closed ladder must obey."""
    for (a, b) in m:
        if abs(b - a) < 4:
            return False
    for e, f in itertools.combinations(m, 2):
        if rel(e, f) == 'D':
            continue
        a, b = min(e), max(e)
        c, d = min(f), max(f)
        if abs(c - a) + abs(d - b) < 3:
            return False
    return True

# ---------- girth, independent of any span/L1 reasoning ----------
def girth(adj):
    n = len(adj)
    best = 10 ** 9
    for s in range(n):
        dist = [-1] * n
        par = [-1] * n
        dist[s] = 0
        Q = deque([s])
        while Q:
            u = Q.popleft()
            if 2 * dist[u] >= best:
                break
            for v in adj[u]:
                if dist[v] < 0:
                    dist[v] = dist[u] + 1
                    par[v] = u
                    Q.append(v)
                elif v != par[u]:
                    c = dist[u] + dist[v] + 1
                    if c < best:
                        best = c
    return best

def ladder_adj(twoN, chords):
    adj = [[] for _ in range(twoN)]
    E = set()
    for i in range(twoN):
        j = (i + 1) % twoN
        E.add(frozenset((i, j)))
    for (a, b) in chords:
        E.add(frozenset((a - 1, b - 1)))
    if len(E) != twoN + len(chords):
        return None            # a chord duplicated a host edge
    for e in E:
        u, v = tuple(e)
        adj[u].append(v)
        adj[v].append(u)
    return adj

# ---------- crossing graph and blocks ----------
def crossing_components(chords):
    k = len(chords)
    adj = [[] for _ in range(k)]
    for i, j in itertools.combinations(range(k), 2):
        if rel(chords[i], chords[j]) == 'X':
            adj[i].append(j)
            adj[j].append(i)
    seen = [False] * k
    comps = []
    for i in range(k):
        if seen[i]:
            continue
        st = [i]
        seen[i] = True
        c = []
        while st:
            u = st.pop()
            c.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    st.append(v)
        comps.append(sorted(c))
    return comps

def blocks(twoN, chords):
    """all proper self-matched circular arcs, as (start, length)"""
    mate = {}
    for (a, b) in chords:
        mate[a] = b
        mate[b] = a
    out = []
    for p in range(1, twoN + 1):
        for L in range(2, twoN, 2):
            arc = set(((p - 1 + t) % twoN) + 1 for t in range(L))
            if all(mate[v] in arc for v in arc):
                out.append((p, L))
    return out

# ---------- 1. how small can a self-matched arc be ----------
print("legal open ladders (= possible self-matched arcs) by size")
W = {}
for k in range(3, 7):
    n = 2 * k
    c = sum(1 for m in matchings(list(range(1, n + 1))) if open_legal(m))
    W[n] = c
    print("   2m = %2d   legal matchings = %d" % (n, c))
print()

# ---------- 2. the equivalence disconnected <=> proper block ----------
print("disconnected crossing graph  <=>  proper self-matched arc")
for n in (6, 8, 10, 12):
    bad = 0
    tot = 0
    for m in matchings(list(range(1, n + 1))):
        tot += 1
        dis = len(crossing_components(m)) > 1
        blk = len(blocks(n, m)) > 0
        if dis != blk:
            bad += 1
    print("   2N = %2d   %6d matchings   mismatches = %d" % (n, tot, bad))
print()

# ---------- 3. the 2N = 24 witness ----------
L1 = [(1, 5), (2, 7), (3, 9), (4, 11), (6, 10), (8, 12)]
L2 = [(a + 12, b + 12) for (a, b) in L1]
CH = L1 + L2
adj = ladder_adj(24, CH)
print("witness at 2N = 24")
print("   chords", CH)
print("   girth (BFS, independent) =", girth(adj))
print("   degrees =", sorted(set(len(a) for a in adj)))
comps = crossing_components(CH)
print("   crossing-graph components =", len(comps))
for c in comps:
    print("      ", [CH[i] for i in c])
bl = blocks(24, CH)
print("   proper self-matched arcs (start, length):", bl)
sav = set()
for (a, b) in CH:
    s = b - a
    sav.add(s - 1)
    sav.add(24 - s - 1)
d = 0
for x in sav:
    d = gcd(d, x)
print("   savings", sorted(sav), " d = gcd =", d)

# each block's own crossing graph
for (p, L) in bl:
    arc = set(((p - 1 + t) % 24) + 1 for t in range(L))
    sub = [e for e in CH if e[0] in arc]
    print("   block at %d len %d: chords %d, own components %d"
          % (p, L, len(sub), len(crossing_components(sub))))

# ---------- 4. cycle spectrum of the witness ----------
def spectrum(adj):
    n = len(adj)
    lens = set()
    for s in range(n):
        # simple cycles whose least vertex is s
        stack = [(s, [s], {s})]
        while stack:
            u, path, used = stack.pop()
            for v in adj[u]:
                if v == s and len(path) > 2:
                    lens.add(len(path))
                elif v > s and v not in used:
                    stack.append((v, path + [v], used | {v}))
    return sorted(lens)

sp = spectrum(adj)
print("   cycle spectrum:", sp)
print("   powers of two present:", [x for x in sp if x & (x - 1) == 0])
