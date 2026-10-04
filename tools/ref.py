"""Independent reference implementation, deliberately using a DIFFERENT
algorithm from egcore.h (subset DP / Held-Karp instead of DFS) so that
agreement is real evidence of correctness.

Also: random cubic graph generation + graph6 encoding.
"""
import random, sys
from itertools import combinations


# ---------- graph6 ----------
def to_graph6(n, edges):
    adj = [[0] * n for _ in range(n)]
    for u, v in edges:
        adj[u][v] = adj[v][u] = 1
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(adj[i][j])
    while len(bits) % 6:
        bits.append(0)
    out = []
    if n <= 62:
        out.append(chr(n + 63))
    else:
        out.append(chr(126))
        out.append(chr(((n >> 12) & 63) + 63))
        out.append(chr(((n >> 6) & 63) + 63))
        out.append(chr((n & 63) + 63))
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return "".join(out)


def from_graph6(s):
    p = 0
    if ord(s[0]) == 126:
        p = 1
        n = 0
        for i in range(3):
            n = (n << 6) | (ord(s[p]) - 63); p += 1
    else:
        n = ord(s[0]) - 63; p = 1
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


# ---------- reference cycle counting via subset DP ----------
def cycle_spectrum_counts(n, edges):
    """Return dict L -> number of simple cycles of length L. Exact.
    Algorithm: for each root r (= minimum vertex of the cycle), DP over
    subsets of {r,...,n-1}: f[S][v] = # simple paths r->v with vertex set S.
    Distinct from the DFS used in the C code."""
    nbr = [0] * n
    for u, v in edges:
        nbr[u] |= 1 << v
        nbr[v] |= 1 << u
    counts = {}
    for r in range(n):
        high = 0
        for v in range(r, n):
            high |= 1 << v
        # f maps (subset, endpoint) -> count; iterate subsets in increasing order
        f = {}
        f[(1 << r, r)] = 1
        # process subsets by popcount order using a dict keyed by subset
        bysize = {}
        bysize.setdefault(1, {})[(1 << r, r)] = 1
        maxsize = n - r
        for size in range(1, maxsize + 1):
            cur = bysize.get(size)
            if not cur:
                continue
            for (S, v), c in cur.items():
                avail = nbr[v] & high & ~S
                m = avail
                while m:
                    w = (m & -m).bit_length() - 1
                    m &= m - 1
                    key = (S | (1 << w), w)
                    d = bysize.setdefault(size + 1, {})
                    d[key] = d.get(key, 0) + c
                # close the cycle
                if size >= 3 and (nbr[v] >> r) & 1:
                    counts[size] = counts.get(size, 0) + c
    # each cycle counted twice (two directions), except length-2 which can't occur
    return {L: c // 2 for L, c in counts.items()}


def pow2_mask(n, edges):
    sp = cycle_spectrum_counts(n, edges)
    m = 0
    k = 2
    while (1 << k) <= n:
        if sp.get(1 << k, 0) > 0:
            m |= 1 << k
        k += 1
    return m, sp


# ---------- random cubic graph (pairing model with rejection) ----------
def random_cubic(n, rng, tries=4000):
    for _ in range(tries):
        stubs = []
        for v in range(n):
            stubs += [v, v, v]
        rng.shuffle(stubs)
        edges = set()
        ok = True
        for i in range(0, len(stubs), 2):
            u, v = stubs[i], stubs[i + 1]
            if u == v or (min(u, v), max(u, v)) in edges:
                ok = False
                break
            edges.add((min(u, v), max(u, v)))
        if ok:
            return sorted(edges)
    return None


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "gen"
    if mode == "gen":
        n = int(sys.argv[2]); cnt = int(sys.argv[3]); seed = int(sys.argv[4])
        rng = random.Random(seed)
        got = 0
        while got < cnt:
            e = random_cubic(n, rng)
            if e is None:
                continue
            print(to_graph6(n, e))
            got += 1
    elif mode == "check":
        # stdin: graph6 lines -> print "<g6> mask=<m> spec=<L:count,...>"
        for line in sys.stdin:
            line = line.strip()
            if not line:
                continue
            n, e = from_graph6(line)
            m, sp = pow2_mask(n, e)
            spec = ",".join("%d:%d" % (L, sp[L]) for L in sorted(sp))
            print("%s mask=%d spec=%s" % (line, m, spec))
