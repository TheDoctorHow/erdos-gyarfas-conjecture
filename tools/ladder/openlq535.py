import sys
sys.setrecursionlimit(10000)

def legal_pair(e, f):
    a, b = e
    c, d = f
    if b < c or d < a:                 # disjoint
        return True
    return abs(c - a) + abs(d - b) >= 3

def gen(n):
    """legal open ladders on [1..n]: spans >= 4, non-disjoint L1 >= 3"""
    mate = [0] * (n + 2)
    cur = []
    out = []

    def rec():
        a = 1
        while a <= n and mate[a]:
            a += 1
        if a > n:
            out.append(tuple(cur))
            return
        for b in range(a + 4, n + 1):
            if mate[b]:
                continue
            e = (a, b)
            if all(legal_pair(e, f) for f in cur):
                mate[a] = b
                mate[b] = a
                cur.append(e)
                rec()
                cur.pop()
                mate[a] = mate[b] = 0
    rec()
    return out

def adj_of(n, ch):
    adj = [[] for _ in range(n + 1)]
    for i in range(1, n):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    for (a, b) in ch:
        adj[a].append(b)
        adj[b].append(a)
    return adj

def has_cycle_len(adj, n, k):
    """is there a simple cycle of length exactly k?"""
    used = [False] * (n + 1)

    def go(start, u, depth):
        if depth == k:
            return start in adj[u]
        for v in adj[u]:
            if v > start and not used[v]:
                used[v] = True
                if go(start, v, depth + 1):
                    used[v] = False
                    return True
                used[v] = False
        return False

    for s in range(1, n + 1):
        used[s] = True
        if go(s, s, 1):
            used[s] = False
            return True
        used[s] = False
    return False

def pow2s(n):
    t = 8
    out = []
    while t <= n:
        out.append(t)
        t *= 2
    return out

for n in (12, 14, 16, 18):
    L = gen(n)
    P = pow2s(n)
    bad = []
    for ch in L:
        adj = adj_of(n, ch)
        if not any(has_cycle_len(adj, n, p) for p in P):
            bad.append(ch)
    print("2m = %2d : %8d legal open ladders, targets %s, WITHOUT a "
          "power-of-two cycle: %d" % (n, len(L), P, len(bad)))
    for ch in bad[:4]:
        print("      ", ch)
    sys.stdout.flush()
