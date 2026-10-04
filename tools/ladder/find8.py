"""List all 8-cycles of cand96.g6 in design coordinates."""
import sys
sys.path.insert(0, ".")
from build96 import edges, N, iv  # noqa: builds graph as side effect

nbr = [[] for _ in range(N)]
for a, b in edges:
    nbr[a].append(b)
    nbr[b].append(a)

def name(v):
    if v < 6:
        return ("u1","u2","u3","v1","v2","v3")[v]
    v -= 6
    c, k = divmod(v, 10)
    x, y = divmod(c, 3)
    return "c%d%d@%d" % (x + 1, y + 1, k + 1)

L = 8
found = set()
def dfs(path, used):
    cur = path[-1]
    k = len(path)
    if k == L:
        if path[0] in nbr[cur]:
            key = frozenset(path)
            if key not in found:
                found.add(key)
                print(" - ".join(name(v) for v in path))
        return
    for w in nbr[cur]:
        if w <= path[0] or used[w]:
            continue
        used[w] = True
        path.append(w)
        dfs(path, used)
        path.pop()
        used[w] = False

used = [False] * N
for r in range(N):
    used[r] = True
    dfs([r], used)
    used[r] = False
print("total 8-cycles:", len(found))
