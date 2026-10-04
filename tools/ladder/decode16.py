"""Decode the 16-cycles of the first stage-1 near-miss in design coordinates."""
import subprocess
from collections import deque

def name(v):
    if v < 6:
        return ("u1","u2","u3","v1","v2","v3")[v]
    v -= 6
    c, k = divmod(v, 10)
    x, y = divmod(c, 3)
    return "c%d%d@%d" % (x + 1, y + 1, k + 1)

def from_g6(s):
    p = 0
    assert ord(s[0]) == 126
    n = 0
    for i in range(1, 4):
        n = (n << 6) | (ord(s[i]) - 63)
    p = 4
    bitpos = 0
    edges = []
    for j in range(1, n):
        for i in range(j):
            byte = bitpos // 6
            bit = 5 - (bitpos % 6)
            if ((ord(s[p + byte]) - 63) >> bit) & 1:
                edges.append((i, j))
            bitpos += 1
    return n, edges

r = subprocess.run(["./egtool.exe", "-q"], stdin=open("b1.g6"),
                   capture_output=True, text=True)
g6 = None
for ln in r.stdout.splitlines():
    if ln.startswith("NEARMISS"):
        g6 = ln.split()[1]
        break
assert g6, "no nearmiss found"
n, edges = from_g6(g6)
nbr = [[] for _ in range(n)]
for a, b in edges:
    nbr[a].append(b)
    nbr[b].append(a)

L = 16
found = []
def bfs_dist(r0):
    d = [10**6] * n
    d[r0] = 0
    q = deque([r0])
    while q:
        u = q.popleft()
        for w in nbr[u]:
            if w >= r0 and d[w] > d[u] + 1:
                d[w] = d[u] + 1
                q.append(w)
    return d

import sys
sys.setrecursionlimit(10000)
def dfs(path, used, dist):
    if len(found) >= 6:
        return
    cur = path[-1]
    k = len(path)
    if k == L:
        if path[0] in nbr[cur] and path[1] < cur:
            found.append(path[:])
        return
    if dist[cur] > L - k:
        return
    for w in nbr[cur]:
        if w <= path[0] or used[w]:
            continue
        used[w] = True
        path.append(w)
        dfs(path, used, dist)
        path.pop()
        used[w] = False

for r0 in range(n):
    if len(found) >= 6:
        break
    dist = bfs_dist(r0)
    used = [False] * n
    used[r0] = True
    dfs([r0], used, dist)

print("sample 16-cycles (up to 6):")
for c in found:
    print("  " + " ".join(name(v) for v in c))
