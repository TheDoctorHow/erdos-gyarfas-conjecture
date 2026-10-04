"""ladder451.py -- THE LADDER BOUND, stress-tested.

ladder450 found ear interiors bounded at ~10 vertices when EVERY
interior vertex carries a self-chord.  Two assumptions must be
relaxed before the bound can be believed:

 (1) UNMATCHED VERTICES.  A vertex whose third edge LEAVES the
     interior (a cross edge) is unmatched.  Allowing x of them can
     only lengthen the ladder.  Measured: max s as a function of x.
     If s <= B (x + 1) then a ladder of length s NEEDS
     x >= s/B - 1 cross edges - a LINEAR cross-edge supply, which
     is exactly what the escape regime was trying to avoid.
 (2) GIRTH.  ladder450 assumed girth >= 5 (span >= 4).  Triangles
     are not themselves forbidden by the conjecture (3 is not a
     power of 2), so the permissive run also allows span 2.
     (Span 3 is always banned: it makes a 4-cycle.)

Forbidden ladder-cycle lengths (unchanged, and justified in 450):
== 0 mod 2^t (progress, or an island member shorter than C) and
every power of 2.  The cycle-length cap only makes the search MORE
permissive, so an upper bound found here is an upper bound on the
truth.
"""
import sys, time
sys.setrecursionlimit(100000)

def forbidden_set(t, cap):
    bad = set()
    M = 1 << t
    for v in range(M, cap + 1, M): bad.add(v)
    p = 2
    while p <= cap: bad.add(p); p *= 2
    return bad

def search(t, mu, smax, xmax, min_span=4, cap=None):
    """Longest prefix 1..s in which every vertex is either matched by
    a self-chord (span in [min_span, mu], span 3 always banned) or is
    one of at most xmax unmatched (cross-edge) vertices."""
    if cap is None: cap = 8 * mu + 12
    bad = forbidden_set(t, cap)
    spans = [z for z in range(min_span, mu + 1) if z != 3]
    best = [0, None, None]
    adj = {}
    for z in range(1, smax + 2):
        adj.setdefault(z, []); adj.setdefault(z + 1, [])
        adj[z].append(z + 1); adj[z + 1].append(z)

    def cycles_through(p, q):
        out = []
        def dfs(v, seen, ln):
            if ln > cap: return
            for w in adj[v]:
                if w == q:
                    if ln + 1 >= 2: out.append(ln + 2)
                elif w not in seen:
                    seen.add(w); dfs(w, seen, ln + 1); seen.discard(w)
        dfs(p, {p}, 0)
        return out

    matched = set(); chords = []; unm = []
    def rec(p, x):
        if p - 1 > best[0]:
            best[0] = p - 1; best[1] = list(chords); best[2] = list(unm)
        if p > smax: return True
        if p in matched: return rec(p + 1, x)
        # option A: match p
        for sp in spans:
            q = p + sp
            if q > smax or q in matched: continue
            if all(ln not in bad for ln in cycles_through(p, q)):
                adj[p].append(q); adj[q].append(p)
                matched.add(p); matched.add(q); chords.append((p, q))
                if rec(p + 1, x): return True
                chords.pop(); matched.discard(p); matched.discard(q)
                adj[p].remove(q); adj[q].remove(p)
        # option B: leave p unmatched (its third edge is a cross edge)
        if x < xmax:
            unm.append(p)
            if rec(p + 1, x + 1): return True
            unm.pop()
        return False

    unbounded = rec(1, 0)
    return best[0], unbounded, best[1], best[2]

print('=== (1) THE LADDER BOUND vs THE CROSS-EDGE BUDGET ===',
      flush=True)
print('girth >= 5 (spans 4..mu).  s = longest interior prefix.')
print()
print('  t  mu    x=0   x=1   x=2   x=3   x=4     growth per extra'
      ' cross edge')
print('-' * 74, flush=True)
for t in (3, 4, 5):
    for mu in (8, 12, 16):
        row = []
        for x in (0, 1, 2, 3, 4):
            s, unb, _, _ = search(t, mu, smax=70, xmax=x)
            row.append(('>=%d' % s) if unb else str(s))
        deltas = []
        vals = []
        for r in row:
            vals.append(int(r.replace('>=', '')))
        for i in range(1, len(vals)): deltas.append(vals[i] - vals[i-1])
        print('  %d  %2d   %5s %5s %5s %5s %5s     %s'
              % (t, mu, row[0], row[1], row[2], row[3], row[4],
                 deltas), flush=True)

print(); print('=== (2) PERMISSIVE RUN: triangles allowed (span >= 2) ==='
               , flush=True)
print('  t  mu    x=0   x=1   x=2   x=3   x=4', flush=True)
for t in (3, 4, 5):
    for mu in (8, 16):
        row = []
        for x in (0, 1, 2, 3, 4):
            s, unb, _, _ = search(t, mu, smax=70, xmax=x, min_span=2)
            row.append(('>=%d' % s) if unb else str(s))
        print('  %d  %2d   %5s %5s %5s %5s %5s'
              % (t, mu, row[0], row[1], row[2], row[3], row[4]),
              flush=True)

print(); print('=== (3) LARGE SPANS: does mu help? ===', flush=True)
for t in (3, 4):
    row = []
    for mu in (16, 24, 32):
        t0 = time.time()
        s, unb, _, _ = search(t, mu, smax=70, xmax=0)
        row.append('mu=%d: %s%d (%.0fs)'
                   % (mu, '>=' if unb else '', s, time.time() - t0))
    print('  t=%d  %s' % (t, '   '.join(row)), flush=True)

print(); print('=== (4) WHAT THE INTERIOR NEEDS ===', flush=True)
for t in (3, 4, 5):
    for m in (3, 5):
        L = m * (1 << t)
        print('  t=%d m=%d: L=%3d, every ear has d > L/2, so the'
              ' interior needs s >= %d' % (t, m, L, L // 2), flush=True)
