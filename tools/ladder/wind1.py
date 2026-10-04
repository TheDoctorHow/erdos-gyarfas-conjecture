"""Three explicit cycles, one per case, covering every girth-legal U(2N,s).

  s = 5                  : S = {0,2,4}          winding 0, length 8
  2N = 2s+4              : one-lap k=2, runs +5,-1              length 8
  2N >= 2s+6, s >= 7     : S = {0,4,s-1,s+5}    winding 0, length 16

Everything is checked by SIMULATION in U(2N,s): walk the cycle, confirm
each step is an edge, the vertices are distinct, and the length is right.
"""
import sys
from wind0 import pairing_is_one_cycle, length_of, cycle_from

def mate_of(n, s):
    mate = [None]*n
    for i in range(0, n, 2):
        mate[i] = (i+s) % n; mate[(i+s) % n] = i
    return mate

def walk_edges(edges, n, mate):
    """Verify: every edge is a real edge of U(n,s), every vertex has degree
    2, and the whole thing is ONE simple cycle.  Returns its length."""
    for u, v in edges:
        if not ((v-u) % n == 1 or (u-v) % n == 1 or mate[u] == v):
            return None
    nb = {}
    for u, v in edges:
        if u == v: return None
        nb.setdefault(u, []).append(v); nb.setdefault(v, []).append(u)
    if any(len(x) != 2 for x in nb.values()): return None
    if len(nb) != len(edges): return None
    start = edges[0][0]; prev = None; cur = start; seen = []
    while True:
        seen.append(cur)
        w = nb[cur][0] if nb[cur][0] != prev else nb[cur][1]
        prev, cur = cur, w
        if cur == start: break
    if len(seen) != len(edges) or len(set(seen)) != len(seen): return None
    return len(seen)

def wind0_edges(S, s, n):
    rungs, host, _ = cycle_from(S, s)
    e = [(m % n, (m+s) % n) for m in rungs]
    e += [(t % n, (t+1) % n) for t in host]
    return e

def onelap_octagon_edges(s, n):
    """k = 2, tails 0 and s+5, runs +5 then -1; used when 2N = 2s+4."""
    seq = [0, s]
    for t in range(1, 6): seq.append((s+t) % n)
    seq.append((2*s+5) % n)          # chord from the tail s+5; 2s+4 = 0 = start
    return [(seq[i], seq[(i+1) % len(seq)]) for i in range(len(seq))]

def construction(s, n):
    if s == 5:
        return 'A', wind0_edges([0, 2, 4], s, n)
    if n == 2*s+4:
        return 'B', onelap_octagon_edges(s, n)
    return 'C', wind0_edges([0, 4, s-1, s+5], s, n)

if __name__ == '__main__':
    smax = int(sys.argv[1]) if len(sys.argv) > 1 else 201
    tot = 0; bad = []; hist = {}
    for s in range(5, smax+1, 2):
        for n in range(2*(s+2), 24*s+1, 2):
            tot += 1
            case, edges = construction(s, n)
            L = walk_edges(edges, n, mate_of(n, s))
            if L is None or L not in (8, 16):
                bad.append((s, n, case, L))
            else:
                hist[(case, L)] = hist.get((case, L), 0) + 1
    print("girth-legal uniform pairs (s odd 5..%d, 2N = 2s+4 .. 24s): %d"
          % (smax, tot))
    print("built and walked in the graph, simple cycle of length 8 or 16: %d"
          % (tot - len(bad)))
    print("failures: %d" % len(bad), bad[:8])
    print("by case and length:", dict(sorted(hist.items())))
