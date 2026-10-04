"""The COMPLETABLE MARGIN.

margin(592) counts the largest girth-legal partial matching with no
L-cycle.  That maximum can be attained by a configuration whose missing
chords are blocked by GIRTH, which says nothing about the target.  The
quantity that measures distance to a counterexample is

  completable margin = N - max { |M| : M girth-legal, no L-cycle, and M
                                 extends to a girth-legal closed ladder }

and it equals

  min over girth-legal closed ladders C of
      (minimum number of chords of C whose deletion destroys every
       L-cycle of C),

because the completable avoiding configurations are exactly the
sub-matchings of closed ladders, and deleting chords only removes cycles.
It is zero exactly when some closed ladder has no L-cycle at all.

Per ladder the value is computed EXACTLY: enumerate every L-cycle by the
(S,eps) dual, then solve the minimum hitting set by branch and bound.
"""
import sys, time, random
from pushL import short, bigpow

def all_L_cycles(n, chords, L):
    """chord-index sets of every cycle of length exactly L"""
    ch = [(min(a, b), max(a, b)) for a, b in chords]
    lower = {a: i for i, (a, b) in enumerate(ch)}
    upper = {b: i for i, (a, b) in enumerate(ch)}
    sel = [False] * len(ch)
    out = []

    def finish(xs):
        edges = [ch[i] for i, s in enumerate(sel) if s]
        idxs = frozenset(i for i, s in enumerate(sel) if s)
        for v in range(n):
            if xs[v]:
                edges.append((v, (v + 1) % n))
        if len(edges) != L:
            return
        adj = {}
        for a, b in edges:
            adj.setdefault(a, []).append(b)
            adj.setdefault(b, []).append(a)
        if any(len(x) != 2 for x in adj.values()):
            return
        st = [next(iter(adj))]; seen = {st[0]}
        while st:
            u = st.pop()
            for w in adj[u]:
                if w not in seen:
                    seen.add(w); st.append(w)
        if len(seen) == len(adj):
            out.append(idxs)

    def scan(eps):
        xs = [0] * n

        def rec(v, t, tot):
            if tot > L:
                return
            if v == n:
                if t == eps:
                    finish(xs)
                return
            if v in lower:
                i = lower[v]
                xs[v] = t
                rec(v + 1, t, tot + t)
                sel[i] = True
                nt = t ^ 1
                xs[v] = nt
                rec(v + 1, nt, tot + nt + 1)
                sel[i] = False
                return
            if v in upper:
                i = upper[v]
                nt = t ^ 1 if sel[i] else t
                xs[v] = nt
                rec(v + 1, nt, tot + nt)
                return
            xs[v] = t
            rec(v + 1, t, tot + t)

        rec(0, eps, 0)

    scan(0); scan(1)
    return out

def min_hitting_set(sets, ub):
    """smallest set of chord indices meeting every set; ub = cut-off"""
    sets = [set(s) for s in sets]
    best = [ub]

    def rec(rem, chosen):
        if not rem:
            best[0] = min(best[0], len(chosen))
            return
        if len(chosen) + 1 > best[0]:
            return
        # branch on the smallest remaining set
        s = min(rem, key=len)
        for c in s:
            rec([t for t in rem if c not in t], chosen | {c})

    rec(sets, frozenset())
    return best[0]

def random_ladder(n, rng, tries=6000):
    for _ in range(tries):
        adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
        mate = [-1] * n
        ok = True
        for v in range(n):
            if mate[v] != -1:
                continue
            cand = [w for w in range(v + 1, n) if mate[w] == -1
                    and not short(adj, n, v, w)]
            if not cand:
                ok = False; break
            w = rng.choice(cand)
            adj[v].append(w); adj[w].append(v)
            mate[v] = w; mate[w] = v
        if ok:
            return [(v, mate[v]) for v in range(n) if mate[v] > v]
    return None

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]))
    trials = int(sys.argv[2])
    for n in [int(x) for x in sys.argv[3:]]:
        L = bigpow(n)
        if L == n:
            print("  2N=%2d trivial" % n); continue
        t = time.time(); best = 10 ** 9; bw = None; ncyc = []
        for _ in range(trials):
            ch = random_ladder(n, rng)
            if ch is None:
                continue
            cyc = all_L_cycles(n, ch, L)
            ncyc.append(len(cyc))
            if not cyc:
                best = 0; bw = ch; break
            h = min_hitting_set(cyc, best)
            if h < best:
                best = h; bw = ch
        print("  2N=%2d  N=%2d  L=%2d   COMPLETABLE MARGIN <= %d"
              "   [%d ladders, L-cycles per ladder %d..%d, %.0fs]"
              % (n, n // 2, L, best, len(ncyc), min(ncyc), max(ncyc),
                 time.time() - t))
        print("        best ladder:", bw)
        sys.stdout.flush()
