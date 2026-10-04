import itertools, sys
from math import gcd
from block534 import matchings, open_legal

def side_graph(m2, ch):
    """host PATH 1..m2 plus the block's chords"""
    adj = [[] for _ in range(m2 + 1)]
    for i in range(1, m2):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    for (a, b) in ch:
        adj[a].append(b)
        adj[b].append(a)
    return adj

def path_lengths(adj, s, t):
    """all lengths of simple s-t paths"""
    out = set()
    n = len(adj) - 1
    used = [False] * (n + 1)

    def go(u, L):
        if u == t:
            out.add(L)
            return
        for v in adj[u]:
            if not used[v]:
                used[v] = True
                go(v, L + 1)
                used[v] = False

    used[s] = True
    go(s, 0)
    return out

def cycle_lengths_side(adj, m2):
    """lengths of cycles living inside one side"""
    out = set()
    for s in range(1, m2 + 1):
        used = [False] * (m2 + 1)
        used[s] = True

        def go(u, L, first):
            for v in adj[u]:
                if v == s and L >= 2 and first != v:
                    out.add(L + 1)
                elif v > s and not used[v]:
                    used[v] = True
                    go(v, L + 1, first)
                    used[v] = False
        for w in adj[s]:
            if w > s:
                used[w] = True
                go(w, 1, w)
                used[w] = False
    return out

def descents(m2, ch):
    adj = side_graph(m2, ch)
    P = path_lengths(adj, 1, m2)
    return set(m2 - 1 - p for p in P), sorted(P)

B = {}
for k in (6, 7):
    m2 = 2 * k
    B[m2] = [m for m in matchings(list(range(1, m2 + 1))) if open_legal(m)]
    print("legal %d-blocks: %d" % (m2, len(B[m2])))
print()

D = {}
for m2 in (12, 14):
    D[m2] = [descents(m2, ch)[0] for ch in B[m2]]
    print("descent sets D(M) for 2m = %d" % m2)
    seen = {}
    for d in D[m2]:
        isint = (sorted(d) == list(range(min(d), max(d) + 1)))
        key = (min(d), max(d), isint, len(set(x % 2 for x in d)))
        seen[key] = seen.get(key, 0) + 1
    for key in sorted(seen):
        lo, hi, isint, npar = key
        print("   min %d  max %2d  interval=%-5s parities=%d   x%d"
              % (lo, hi, isint, npar, seen[key]))
    print("   every D(M) = [0, max] :",
          all(sorted(d) == list(range(0, max(d) + 1)) for d in D[m2]))
    print("   min over blocks of max D:", min(max(d) for d in D[m2]))
    print()

def covers(D1, D2, target):
    return any(target - x in D2 for x in D1)

print("2N = 24, disconnected: need 8 in D(M1) + D(M2)")
bad = [(i, j) for i in range(len(D[12])) for j in range(len(D[12]))
       if not covers(D[12][i], D[12][j], 8)]
print("   failing pairs:", len(bad), "of", len(D[12]) ** 2)

print("2N = 26, disconnected: need 10 in D(M12) + D(M14)")
bad2 = [(i, j) for i in range(len(D[12])) for j in range(len(D[14]))
        if not covers(D[12][i], D[14][j], 10)]
print("   failing pairs:", len(bad2), "of", len(D[12]) * len(D[14]))
print()

c8 = sum(1 for ch in B[12]
         if 8 in cycle_lengths_side(side_graph(12, ch), 12))
print("12-blocks whose own side already holds an 8-cycle: %d of %d"
      % (c8, len(B[12])))
