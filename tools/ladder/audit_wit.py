"""Audit a partial-matching witness:
  - girth (BFS, independent of `short`)
  - no cycle of length L (dual detector of 581)
  - WHY it is maximal: for every unmatched pair, does adding that chord
    break girth, or create an L-cycle?
  - is it COMPLETABLE, i.e. does the unmatched set admit a girth-legal
    perfect matching?  A configuration that is only saturated because its
    leftover vertices are too close together says nothing about the
    target.
"""
import sys
from dualdet import detect

def girth(n, chords):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    for a, b in chords:
        adj[a].append(b); adj[b].append(a)
    best = 10 ** 9
    for s in range(n):
        dist = [-1] * n; par = [-1] * n
        dist[s] = 0; q = [s]; i = 0
        while i < len(q):
            u = q[i]; i += 1
            for w in adj[u]:
                if dist[w] < 0:
                    dist[w] = dist[u] + 1; par[w] = u; q.append(w)
                elif w != par[u]:
                    best = min(best, dist[u] + dist[w] + 1)
    return best

def legal_add(n, chords, a, b):
    """would adding (a,b) keep girth >= 5?  distance >= 4 in the graph"""
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    for x, y in chords:
        adj[x].append(y); adj[y].append(x)
    dist = {a: 0}; q = [a]; i = 0
    while i < len(q):
        u = q[i]; i += 1
        if dist[u] >= 3:
            continue
        for w in adj[u]:
            if w not in dist:
                dist[w] = dist[u] + 1; q.append(w)
    return dist.get(b, 99) >= 4

def completable(n, chords):
    """does the unmatched set admit a girth-legal perfect matching?"""
    cov = set()
    for a, b in chords:
        cov.add(a); cov.add(b)
    free = sorted(v for v in range(n) if v not in cov)
    cur = list(chords)

    def rec(rem):
        if not rem:
            return list(cur)
        a = rem[0]
        for j in range(1, len(rem)):
            b = rem[j]
            if not legal_add(n, cur, a, b):
                continue
            cur.append((a, b))
            r = rec([x for x in rem[1:] if x != b])
            if r:
                return r
            cur.pop()
        return None

    return rec(free)

if __name__ == '__main__':
    n = int(sys.argv[1]); L = int(sys.argv[2])
    ch = [(min(a, b), max(a, b)) for a, b in eval(sys.argv[3])]
    cov = set()
    for a, b in ch:
        cov.add(a); cov.add(b)
    free = sorted(v for v in range(n) if v not in cov)
    print("  n=%d  L=%d  chords=%d  unmatched=%s" % (n, L, len(ch), free))
    print("  girth =", girth(n, ch))
    print("  has an L-cycle:", bool(detect(n, ch, L)))
    print("  why maximal -- each addable pair:")
    any_legal = False
    for i in range(len(free)):
        for j in range(i + 1, len(free)):
            a, b = free[i], free[j]
            ok = legal_add(n, ch, a, b)
            if not ok:
                print("      (%d,%d)  BLOCKED BY GIRTH" % (a, b))
            else:
                any_legal = True
                hit = detect(n, ch + [(a, b)], L)
                print("      (%d,%d)  girth ok, creates an L-cycle: %s"
                      % (a, b, bool(hit)))
    c = completable(n, ch)
    print("  completable to a girth-legal closed ladder:", bool(c))
    if c:
        print("      completion:", sorted(c))
        print("      that ladder has an L-cycle:", bool(detect(n, c, L)))
