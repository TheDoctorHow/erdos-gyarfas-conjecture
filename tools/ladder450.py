"""ladder450.py -- THE ESCAPE REGIME: can a long ear interior exist?

The regime (448.6): every ear shift is large, so self-chords have
room and the Disjoint-Chord Kill does not fire.  Then the interior
of an ear is a LADDER: a path v_1..v_s with (nearly) every vertex
carrying one self-chord of span <= mu.

THE SQUEEZE.  C is the SHORTEST island member, so every cycle of G
shorter than L has v_2 < t (a shorter cycle with v_2 = t would beat
C; a cycle with v_2 > t is progress).  Every ladder cycle has
length <= (chords used) x mu + O(1), and in the escape regime
mu ~ sqrt(L) << L.  Hence EVERY LADDER CYCLE LENGTH must satisfy

    length =/= 0 (mod 2^t)      and      length =/= 2^k.

Also, two structural facts already proved:
    - girth >= 5: every chord span sigma >= 4;
    - 431 local minimality on BOTH arcs: d > L/2 for every ear, so
      the interior length is s >= L/2 - hence the ladder must be
      LONG (>= 2^{t-1} m).

QUESTION (decisive for the escape): for given t and mu, how long
can such a ladder be?  If the maximum length is bounded below
2^{t-1} m, the escape regime is EMPTY and J7's last hiding place
closes.

Machine: exhaustive DFS over matchings (chords placed left to
right, spans in [4, mu]), with every cycle through each new chord
enumerated and tested.  Reports the maximum attainable s, and
whether the pattern is periodic (i.e. extends forever).
"""
import sys, time
sys.setrecursionlimit(100000)

def forbidden_set(t, cap):
    """Lengths a ladder cycle may not have."""
    bad = set()
    M = 1 << t
    for v in range(M, cap + 1, M): bad.add(v)     # == 0 mod 2^t
    p = 2
    while p <= cap: bad.add(p); p *= 2            # powers of two
    return bad

def cycles_through(adj, p, q, cap):
    """All simple p-q path lengths in adj (chord (p,q) excluded)."""
    out = []
    def dfs(v, seen, ln):
        if ln > cap: return
        for w in adj[v]:
            if w == q:
                if ln + 1 >= 2: out.append(ln + 1)
            elif w not in seen:
                seen.add(w); dfs(w, seen, ln + 1); seen.discard(w)
    dfs(p, {p}, 0)
    return out

def search(t, mu, smax, cap=None, want_all=False):
    """Longest ladder: path 1..s, matching with spans in [4, mu],
    every vertex matched except possibly the last few.  Returns
    (best_s, witness_chords)."""
    if cap is None: cap = 6 * mu + 8
    bad = forbidden_set(t, cap)
    best = [0, None]
    adj = {}
    def add_edge(a, b):
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    def del_edge(a, b):
        adj[a].remove(b); adj[b].remove(a)

    def ok_new_chord(p, q):
        # all cycles through the new chord (p,q)
        for ln in cycles_through(adj, p, q, cap):
            if (ln + 1) in bad: return False
        return True

    chords = []
    matched = set()

    def rec(p):
        """positions 1..p-1 all decided; extend."""
        # record progress: the ladder is valid up to p-1
        if p - 1 > best[0]:
            best[0] = p - 1; best[1] = list(chords)
        if p > smax: return True          # reached the cap: unbounded
        if p in matched:
            return rec(p + 1)
        # p must be matched to some q in [p+4, p+mu]
        for q in range(p + 4, min(p + mu, smax) + 1):
            if q in matched: continue
            # ensure the path edges up to q exist
            for z in range(p, q):
                if not (z + 1 in adj.get(z, [])):
                    add_edge(z, z + 1)
            if ok_new_chord(p, q):
                add_edge(p, q)
                matched.add(p); matched.add(q); chords.append((p, q))
                if rec(p + 1): return True
                chords.pop(); matched.discard(p); matched.discard(q)
                del_edge(p, q)
        return False

    # seed the path
    for z in range(1, smax + 1):
        add_edge(z, z + 1) if False else None
    # (edges are added lazily inside rec)
    add_edge(1, 2)
    for z in range(2, smax + 1): add_edge(z, z + 1)
    unbounded = rec(1)
    return best[0], best[1], unbounded

print('=== THE LADDER SEARCH ===', flush=True)
print('forbidden cycle lengths: == 0 mod 2^t, and all powers of 2')
print()
print(' t  mu   max ladder length s      needed s = 2^(t-1) m (m=3)'
      '   verdict')
print('-' * 76, flush=True)
for t in (3, 4):
    need = (1 << (t - 1)) * 3
    for mu in (6, 8, 10, 12, 16):
        t0 = time.time()
        s, wit, unb = search(t, mu, smax=max(60, need + 8))
        verdict = ('reaches the cap (extends indefinitely)' if unb
                   else 'BOUNDED at %d' % s)
        print(' %d  %2d   %-24s %-24d %s  [%.0fs]'
              % (t, mu, ('>= %d' % s) if unb else str(s), need,
                 verdict, time.time() - t0), flush=True)
        if not unb and s < need:
            print('      => ESCAPE DEAD at t=%d, mu=%d: the interior'
                  ' cannot reach length %d' % (t, mu, need), flush=True)
        if unb and wit is not None:
            print('      witness prefix chords: %s'
                  % wit[:10], flush=True)
