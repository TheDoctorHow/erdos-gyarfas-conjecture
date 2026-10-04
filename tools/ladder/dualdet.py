"""Independent detector for a cycle of length exactly L in host + chords,
via the (S, eps) dual of 581 -- no path search.

Every element of the cycle space is (S, eps): choose a chord subset S and
a starting parity eps; the host-edge indicator satisfies
    x_v = x_{v-1}  XOR  [v is an endpoint of S],
so the element is forced once S and x_{-1} = eps are fixed, and its
length is |S| + sum_v x_v.  Scanning v = 0..n-1 and branching only at the
LOWER endpoint of each chord gives a search whose partial length is
monotone, so a partial length above L kills the branch.

The element is a simple cycle iff its edge set is connected (every vertex
it touches already has degree exactly 2 by construction).
"""
import sys

def detect(n, chords, L, want_witness=False):
    ch = [(min(a, b), max(a, b)) for a, b in chords]
    lower = {}
    upper = {}
    for i, (a, b) in enumerate(ch):
        lower[a] = i
        upper[b] = i
    sel = [False] * len(ch)
    found = [None]

    def finish(eps, xs):
        edges = []
        for i, s in enumerate(sel):
            if s:
                edges.append(ch[i])
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
        st = [next(iter(adj))]
        seen = set(st)
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
            if found[0] is not None:
                return True
            if tot > L:
                return False
            if v == n:
                if t != eps:
                    return False
                return finish(eps, xs)
            if v in lower:
                i = lower[v]
                # not selected
                xs[v] = t
                if rec(v + 1, t, tot + t):
                    return True
                # selected
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

    for eps in (0, 1):
        if scan(eps):
            return found[0] if want_witness else True
    return None if want_witness else False

if __name__ == '__main__':
    n = int(sys.argv[1]); L = int(sys.argv[2])
    chords = eval(sys.argv[3])
    w = detect(n, chords, L, want_witness=True)
    print("  n=%d  L=%d  chords=%d" % (n, L, len(chords)))
    if w is None:
        print("  NO cycle of length exactly %d" % L)
    else:
        print("  cycle of length %d found:" % L, w)
