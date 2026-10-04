"""Explicit construction refuting the crux lemma of the truncation route.

CLAIM UNDER TEST:  every cubic graph of girth 5 has a cycle of length 6, 7 or 8.

CONSTRUCTION.  Let H be a cubic graph of girth >= 9.  Pick five edges of H that are
pairwise far apart, subdivide each once (creating s_1..s_5 of degree 2), take a new
5-cycle v_1..v_5, and join v_j to s_j.  Every vertex now has degree 3.

  * the new 5-cycle has length 5, so girth <= 5;
  * every cycle lying inside the subdivided H has length >= girth(H) >= 9;
  * every cycle using the pentagon enters at some v_i and leaves at some v_j, so has
    length  1 + d(s_i,s_j) + 1 + a  with a in {1,2}; choosing the five edges pairwise
    far apart makes this >= 9.
  So girth(G) = 5 and G has no cycle of length 6, 7 or 8.
"""
import sys, collections, itertools
sys.path.insert(0, '.')
from ref import to_graph6


def read_conder(path, want_girth, nmax):
    """Yield (name, n, edges) for cubic census graphs with girth >= want_girth."""
    with open(path, 'r') as f:
        for line in f:
            p = line.split()
            if len(p) < 4:
                continue
            name = p[0]
            try:
                n = int(p[1]); m = int(p[3])
            except ValueError:
                continue
            if n > nmax:
                continue
            nums = list(map(int, p[4:4 + 2 * m]))
            if len(nums) != 2 * m or m != 3 * n // 2:
                continue
            edges = [(nums[2 * i], nums[2 * i + 1]) for i in range(m)]
            if girth(n, edges) >= want_girth:
                yield name, n, edges


def adjacency(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v); adj[v].append(u)
    return adj


def girth(n, edges):
    adj = adjacency(n, edges)
    best = 10 ** 9
    for s in range(n):
        dist = {s: 0}; par = {s: -1}
        q = collections.deque([s])
        while q:
            u = q.popleft()
            if dist[u] * 2 >= best:
                break
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1; par[w] = u; q.append(w)
                elif w != par[u]:
                    best = min(best, dist[u] + dist[w] + 1)
    return best


def bfs(adj, s):
    d = {s: 0}
    q = collections.deque([s])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in d:
                d[w] = d[u] + 1; q.append(w)
    return d


def pick_far_edges(n, edges, k, sep):
    """Greedily pick k edges pairwise at vertex-distance >= sep."""
    adj = adjacency(n, edges)
    dist = {v: bfs(adj, v) for v in range(n)}

    def edist(e, f):
        return min(dist[a][b] for a in e for b in f)

    chosen = []
    for e in edges:
        if all(edist(e, f) >= sep for f in chosen):
            chosen.append(e)
            if len(chosen) == k:
                return chosen
    return None


def build(n, edges, five):
    """Subdivide the five chosen edges, attach a pentagon."""
    rest = [e for e in edges if e not in five and (e[1], e[0]) not in five]
    new = list(rest)
    idx = n
    subs = []
    for (u, v) in five:
        s = idx; idx += 1
        new += [(u, s), (s, v)]
        subs.append(s)
    pent = list(range(idx, idx + 5)); idx += 5
    for i in range(5):
        new.append((pent[i], pent[(i + 1) % 5]))
        new.append((pent[i], subs[i]))
    return idx, new


if __name__ == '__main__':
    census = sys.argv[1] if len(sys.argv) > 1 else 'conder_census.txt'
    want = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 130
    sep = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    for name, n, edges in read_conder(census, want, nmax):
        g = girth(n, edges)
        five = pick_far_edges(n, edges, 5, sep)
        if five is None:
            print('# %s n=%d girth=%d : no 5 edges pairwise >= %d apart' % (name, n, g, sep),
                  file=sys.stderr)
            continue
        N, E = build(n, edges, five)
        deg = collections.Counter()
        for u, v in E:
            deg[u] += 1; deg[v] += 1
        assert set(deg.values()) == {3}, deg
        assert len(set(map(lambda e: tuple(sorted(e)), E))) == len(E)
        print(to_graph6(N, [tuple(sorted(e)) for e in E]))
        print('# base %s n=%d girth=%d -> constructed n=%d' % (name, n, g, N), file=sys.stderr)
        break
