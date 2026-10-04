"""dualdet with an admissible suffix bound, so it reaches 2N = 60+.

While the scan is INSIDE (parity t = 1) every host edge is used until the
next toggle, and a toggle can only happen at a chord endpoint.  So from
position v with t = 1 the element must still collect at least
    nextend(v) - v
host edges, where nextend(v) is the least chord endpoint >= v (or n).
Adding that to the running total is admissible, and it is what makes the
length-32 scan finish on a 62-vertex ladder.
"""
import sys, time

def detect2(n, chords, L, want_witness=False, deadline=None,
            eps_list=(0, 1)):
    ch = [(min(a, b), max(a, b)) for a, b in chords]
    lower = {a: i for i, (a, b) in enumerate(ch)}
    upper = {b: i for i, (a, b) in enumerate(ch)}
    ends = set(lower) | set(upper)
    nextend = [n] * (n + 1)
    for v in range(n - 1, -1, -1):
        nextend[v] = v if v in ends else nextend[v + 1]
    sel = [False] * len(ch)
    found = [None]
    nodes = [0]

    def finish(xs):
        edges = [ch[i] for i, s in enumerate(sel) if s]
        for v in range(n):
            if xs[v]:
                edges.append((v, (v + 1) % n))
        if len(edges) != L:
            return False
        adj = {}
        for a, b in edges:
            adj.setdefault(a, []).append(b)
            adj.setdefault(b, []).append(a)
        if any(len(x) != 2 for x in adj.values()):
            return False
        st = [next(iter(adj))]; seen = {st[0]}
        while st:
            u = st.pop()
            for w in adj[u]:
                if w not in seen:
                    seen.add(w); st.append(w)
        if len(seen) != len(adj):
            return False
        found[0] = sorted(edges)
        return True

    def scan(eps):
        xs = [0] * n

        def rec(v, t, tot):
            nodes[0] += 1
            if found[0] is not None:
                return True
            if deadline and not nodes[0] % 200000 and time.time() > deadline:
                raise TimeoutError
            if tot > L:
                return False
            if t and tot + (nextend[v] - v) > L:
                return False
            if v == n:
                return t == eps and finish(xs)
            if v in lower:
                i = lower[v]
                xs[v] = t
                if rec(v + 1, t, tot + t):
                    return True
                sel[i] = True
                nt = t ^ 1
                xs[v] = nt
                if rec(v + 1, nt, tot + nt + 1):
                    return True
                sel[i] = False
                return False
            if v in upper:
                i = upper[v]
                nt = t ^ 1 if sel[i] else t
                xs[v] = nt
                return rec(v + 1, nt, tot + nt)
            xs[v] = t
            return rec(v + 1, t, tot + t)

        return rec(0, eps, 0)

    for eps in eps_list:
        if scan(eps):
            return (found[0] if want_witness else True), nodes[0]
    return (None if want_witness else False), nodes[0]

if __name__ == '__main__':
    sys.setrecursionlimit(200000)
    n = int(sys.argv[1]); L = int(sys.argv[2])
    chords = eval(sys.argv[3])
    t = time.time()
    r, nd = detect2(n, chords, L, want_witness=True)
    print("  n=%d L=%d chords=%d  nodes=%d  [%.0fs]"
          % (n, L, len(chords), nd, time.time() - t))
    print("  " + ("NO cycle of length exactly %d" % L if r is None
                  else "cycle of length %d: %s" % (L, r)))
