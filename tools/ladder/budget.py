"""How many partners does each forbidden length exclude?

Adding a chord (a,b) creates a cycle of length L exactly when a simple
path of length L-1 already joins a to b.  So the vertices excluded by
"no L-cycle" are the endpoints of length-(L-1) paths from a, and there
are at most 3 * 2^(L-2) of them - A BOUND INDEPENDENT OF N.
"""
import sys, random, time

def endpoints(adj, n, a, depths):
    """set of b joined to a by a simple path of length d, for d in depths."""
    dmax = max(depths)
    out = set()
    stack = [(a, 1 << a, 0)]
    while stack:
        u, vis, d = stack.pop()
        if d in depths and u != a: out.add(u)
        if d == dmax: continue
        for w in adj[u]:
            if (vis >> w) & 1: continue
            stack.append((w, vis | (1 << w), d+1))
    return out

def build(n, rng, forbid_lengths, tries=1, verbose=False):
    """greedy closed ladder avoiding cycles of the given lengths (and of
    length <= 4, i.e. girth-legal)."""
    depths = set(L-1 for L in forbid_lengths) | {1, 2, 3}
    for _ in range(tries):
        adj = [[(v+1) % n, (v-1) % n] for v in range(n)]
        mate = [-1]*n
        order = list(range(n)); ok = True
        stats = []
        for a in order:
            if mate[a] != -1: continue
            bad = endpoints(adj, n, a, depths)
            cand = [b for b in range(n)
                    if b != a and mate[b] == -1 and b not in bad]
            stats.append((len(bad), len(cand)))
            if not cand: ok = False; break
            b = rng.choice(cand)
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
        if ok and all(m != -1 for m in mate):
            return mate, stats
        if verbose:
            print("      failed after %d chords; last exclusion set %d"
                  % (len(stats), stats[-1][0] if stats else 0))
    return None, stats

if __name__ == '__main__':
    rng = random.Random(31337)
    print("== the exclusion budget ==")
    print("  forbidding cycles of length L excludes at most 3*2^(L-2) partners:")
    for L in (5, 8, 16, 32):
        print("     L = %2d  ->  <= %d" % (L, 3*2**(L-2)))
    print()
    print("== can a greedy avoid octagons?  and octagons AND hexadecagons? ==")
    for n in (28, 40, 60, 100, 150, 200):
        t = time.time()
        m1, s1 = build(n, rng, (8,), tries=20)
        e1 = max(x[0] for x in s1) if s1 else 0
        t1 = time.time()-t
        t = time.time()
        m2, s2 = build(n, rng, (8, 16), tries=3)
        e2 = max(x[0] for x in s2) if s2 else 0
        print("  2N=%3d :  no-8  %-7s (max excluded %4d, %.0fs)   "
              "no-8-no-16  %-7s (max excluded %5d, %.0fs)"
              % (n, "OK" if m1 else "failed", e1, t1,
                 "OK" if m2 else "failed", e2, time.time()-t))
        sys.stdout.flush()
