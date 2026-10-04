"""Upper bounds on the MARGIN by chord deletion.

margin = N - (max chords in a girth-legal partial matching with no cycle
of length L).  Deleting chords only destroys cycles, so if r chords can
be deleted from a girth-legal CLOSED ladder so that no L-cycle survives,
then max >= N - r and hence margin <= r.

This searches for a small deletion set directly - a hitting set for the
L-cycles - which is far stronger than growing random partial matchings,
because it starts N chords deep instead of climbing to that depth.
"""
import sys, time, random
from collections import deque
from pushL import short, bigpow

def path_through(adj, n, a, b, L):
    """vertex list of a simple L-cycle using edge (a,b), or None."""
    target = L - 1
    d = [-1] * n; d[b] = 0; q = deque([b])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if d[w] < 0:
                d[w] = d[u] + 1; q.append(w)
    if d[a] < 0:
        return None
    seen = bytearray(n); seen[a] = 1
    path = [a]
    def rec(u, dep):
        rem = target - dep
        if d[u] > rem:
            return False
        if rem == 0:
            return u == b
        for w in adj[u]:
            if w == b:
                if rem == 1:
                    path.append(b); return True
                continue
            if seen[w]:
                continue
            seen[w] = 1; path.append(w)
            if rec(w, dep + 1):
                return True
            path.pop(); seen[w] = 0
        return False
    for w in adj[a]:
        if w == b or seen[w]:
            continue
        seen[w] = 1; path.append(w)
        if rec(w, 1):
            return path
        path.pop(); seen[w] = 0
    return None

def find_L_cycle(adj, n, chords, L, rng=None):
    """any L-cycle; every L-cycle (L < n) must use a chord."""
    order = list(range(len(chords)))
    if rng: rng.shuffle(order)
    for i in order:
        a, b = chords[i]
        p = path_through(adj, n, a, b, L)
        if p:
            return p
    return None

def chords_on(path, chordset):
    out = []
    for i in range(len(path)):
        e = (min(path[i], path[i - 1]), max(path[i], path[i - 1]))
        if e in chordset:
            out.append(e)
    return out

def random_ladder(n, rng, g=5, tries=4000):
    for _ in range(tries):
        adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
        mate = [-1] * n
        ok = True
        for v in range(n):
            if mate[v] != -1:
                continue
            cand = [w for w in range(v + 1, n) if mate[w] == -1
                    and not short(adj, n, v, w, g - 1)]
            if not cand:
                ok = False; break
            w = rng.choice(cand)
            adj[v].append(w); adj[w].append(v)
            mate[v] = w; mate[w] = v
        if ok:
            return [(v, mate[v]) for v in range(n) if mate[v] > v]
    return None

def delete_to_avoid(n, L, rng, restarts, verbose=False):
    best = None; bestcfg = None
    for _ in range(restarts):
        ch = random_ladder(n, rng)
        if ch is None:
            continue
        chords = list(ch)
        adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
        for a, b in chords:
            adj[a].append(b); adj[b].append(a)
        cs = set(chords)
        removed = 0
        while True:
            p = find_L_cycle(adj, n, list(cs), L, rng)
            if p is None:
                break
            cand = chords_on(p, cs)
            if not cand:
                break
            a, b = rng.choice(cand)
            adj[a].remove(b); adj[b].remove(a)
            cs.discard((a, b))
            removed += 1
        if best is None or removed < best:
            best = removed; bestcfg = sorted(cs)
            if verbose:
                print("      %d deletions" % removed); sys.stdout.flush()
    return best, bestcfg

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    rng = random.Random(int(sys.argv[1]))
    restarts = int(sys.argv[2])
    for n in [int(x) for x in sys.argv[3:]]:
        L = bigpow(n)
        if L == n:
            print("  2N=%2d trivial" % n); continue
        t = time.time()
        r, cfg = delete_to_avoid(n, L, rng, restarts)
        print("  2N=%2d  N=%2d  L=%2d   deletions = %d   =>  MARGIN <= %d"
              "   [%.0fs]" % (n, n // 2, L, r, r, time.time() - t))
        print("        surviving chords:", cfg)
        sys.stdout.flush()
