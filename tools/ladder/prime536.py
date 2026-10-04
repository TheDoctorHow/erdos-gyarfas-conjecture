import sys
sys.setrecursionlimit(20000)

W36 = [(1,7),(2,11),(3,15),(4,26),(5,32),(6,20),(8,30),(9,18),(10,25),
       (12,34),(13,19),(14,28),(16,31),(17,23),(21,27),(22,36),(24,33),
       (29,35)]

def open_adj(n, ch):
    adj = [[] for _ in range(n + 1)]
    for i in range(1, n):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    for (a, b) in ch:
        adj[a].append(b)
        adj[b].append(a)
    return adj

def count_len(adj, n, k):
    used = [False] * (n + 1)
    tot = 0
    def go(s, u, d):
        nonlocal tot
        if d == k:
            if s in adj[u]:
                tot += 1
            return
        for v in adj[u]:
            if v > s and not used[v]:
                used[v] = True
                go(s, v, d + 1)
                used[v] = False
    for s in range(1, n + 1):
        used[s] = True
        go(s, s, 1)
        used[s] = False
    return tot // 2

def legal_open(n, ch):
    for (a, b) in ch:
        if b - a < 4:
            return False
    for i in range(len(ch)):
        for j in range(i + 1, len(ch)):
            a, b = ch[i]
            c, d = ch[j]
            if b < c or d < a:
                continue
            if abs(c - a) + abs(d - b) < 3:
                return False
    return True

def crossing_conn(ch):
    k = len(ch)
    adj = [[] for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            a, b = ch[i]
            c, d = ch[j]
            if a < c < b < d or c < a < d < b:
                adj[i].append(j)
                adj[j].append(i)
    seen = [False] * k
    st = [0]
    seen[0] = True
    n = 1
    while st:
        u = st.pop()
        for v in adj[u]:
            if not seen[v]:
                seen[v] = True
                n += 1
                st.append(v)
    return n == k

print("THE ARCHIVE'S OWN HARD CASE: 523.0's witness, cut open")
ch = sorted(W36)
adj = open_adj(36, ch)
print("   legal open ladder:", legal_open(36, ch))
print("   PRIME (crossing graph connected):", crossing_conn(ch))
for k in (5, 6, 7, 8, 16):
    print("   %2d-cycles: %d" % (k, count_len(adj, 36, k)))
sys.stdout.flush()
print()

# ---- smallest octagon-free legal open ladder, by pruned search ----
def path_exists(adj, a, b, L, banned):
    """is there a simple path a->b of length exactly L, not using
    the banned edge?"""
    seen = [False] * (len(adj))
    seen[a] = True
    def go(u, d):
        if d == L:
            return u == b
        if u == b:
            return False
        for v in adj[u]:
            if seen[v]:
                continue
            if (u, v) == banned or (v, u) == banned:
                continue
            seen[v] = True
            if go(v, d + 1):
                seen[v] = False
                return True
            seen[v] = False
        return False
    return go(a, 0)

def search_free(n, targets, want_one=True):
    """legal open ladders on [1..n] with no cycle of any length in
    targets.  Prune: a new cycle must use the newly added chord."""
    mate = [0] * (n + 2)
    adj = [[] for _ in range(n + 1)]
    for i in range(1, n):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    cur = []
    found = []

    def ok_pair(e, f):
        a, b = e
        c, d = f
        if b < c or d < a:
            return True
        return abs(c - a) + abs(d - b) >= 3

    def rec():
        if found and want_one:
            return
        a = 1
        while a <= n and mate[a]:
            a += 1
        if a > n:
            found.append(tuple(cur))
            return
        for b in range(a + 4, n + 1):
            if mate[b]:
                continue
            e = (a, b)
            if not all(ok_pair(e, f) for f in cur):
                continue
            adj[a].append(b)
            adj[b].append(a)
            bad = any(path_exists(adj, a, b, t - 1, (a, b)) for t in targets)
            if not bad:
                mate[a] = b
                mate[b] = a
                cur.append(e)
                rec()
                cur.pop()
                mate[a] = mate[b] = 0
            adj[a].pop()
            adj[b].pop()
            if found and want_one:
                return
    rec()
    return found

print("smallest legal open ladder with NO OCTAGON")
for n in range(12, 33, 2):
    f = search_free(n, (8,))
    print("   2m = %2d : %s" % (n, "FOUND " + str(f[0]) if f else "none"))
    sys.stdout.flush()
    if f:
        break
