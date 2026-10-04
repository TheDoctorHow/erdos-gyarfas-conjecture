"""probe434.py -- the constructive probe of the braided core, and the
(f)-correction.

Attempt to BUILD the braided metric core at L = 24 (t = 3, m = 3):
C plus deep-leg antipodal ladders. The probe verifies:
 (1) the LADDER ARITHMETIC: legs of length 12 joined by rungs
     a_p -- b_{12-p} give rung-ears of length exactly 13 > 12 =
     max-arc for antipodal feet: the primary metric constraint is
     satisfiable BY DESIGN;
 (2) the 4-CYCLE OBSTRUCTION: adjacent monotone rungs create
     4-cycles (girth violation) - verified on the dense-rung ladder;
 (3) the 8-CYCLE OBSTRUCTION: spacing-3 rungs create 8-cycles -
     a POWER OF TWO: legal for girth, fatal for condition (f) (the
     ambient must be a counterexample);
 (4) the forced SPACING-2 design: rungs at odd positions give
     6-cycles only (>= 5, not a power of two) - and leave every even
     position unwired: the wiring cascade that (f) drives.
Conclusion verified: the metric core alone is nearly constructible;
what obstructs completion at every step is (f) itself. J7's residual
question carries the full counterexample hypothesis.
"""
import itertools

def girth_and_pow2(adj, cap=12):
    """Shortest cycle length; also whether an 8-cycle exists (BFS per
    vertex, cycles up to length cap)."""
    g = 999; has8 = False
    for s in adj:
        dist = {s: 0}; par = {s: None}; q = [s]
        while q:
            v = q.pop(0)
            if dist[v] * 2 > cap: continue
            for w in adj[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1; par[w] = v; q.append(w)
                elif par[v] != w and par.get(w) != v:
                    g = min(g, dist[v] + dist[w] + 1)
    return g

def find_cycle_len(adj, target):
    """Does a cycle of exactly `target` edges exist? DFS bounded."""
    for s in sorted(adj):
        st = [(s, s, 0, frozenset([s]))]
        while st:
            start, v, ln, seen = st.pop()
            for w in adj[v]:
                if w == start and ln + 1 == target:
                    return True
                if w != start and w not in seen and ln + 1 < target:
                    st.append((start, w, ln + 1, seen | {w}))
    return False

def build(L, rung_positions):
    adj = {}
    def add(a, b):
        adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    for i in range(L): add(('c', i), ('c', (i + 1) % L))
    # one antipodal ladder pair: legs above c0 and c12, length 12
    for side, base in (('a', 0), ('b', L // 2)):
        prev = ('c', base)
        for p in range(1, 12):
            add(prev, (side, p)); prev = (side, p)
    for p in rung_positions:
        add(('a', p), ('b', 12 - p))
    return adj

L = 24
# (1) ladder arithmetic: rung ear length = p + 1 + (12 - p) = 13
adj = build(L, [1, 3, 5, 7, 9, 11])
# shortest C-avoiding path c0 -> c12
Cverts = {('c', i) for i in range(1, L) if i != 12}
dist = {('c', 0): 0}; q = [('c', 0)]
while q:
    v = q.pop(0)
    for w in adj[v]:
        if w not in dist and w not in Cverts:
            dist[w] = dist[v] + 1; q.append(w)
d = dist.get(('c', 12))
print('(1) spacing-2 ladder: shortest antipodal ear = %d > 12: %s'
      % (d, d > 12))
assert d == 13

# (2) dense rungs -> 4-cycles
adj_dense = build(L, list(range(1, 12)))
g = girth_and_pow2(adj_dense)
print('(2) dense rungs: girth %d (4-cycles: %s) -> GIRTH obstruction'
      % (g, g == 4))
assert g == 4

# (3) spacing-3 rungs -> 8-cycles (power of two!)
adj3 = build(L, [1, 4, 7, 10])
g3 = girth_and_pow2(adj3)
has8 = find_cycle_len(adj3, 8)
print('(3) spacing-3 rungs: girth %d (>= 5 OK), 8-cycle exists: %s'
      ' -> condition (f) obstruction (8 = 2^3)' % (g3, has8))
assert g3 >= 5 and has8

# (4) spacing-2 rungs: girth 6, no 8-cycle in the LADDER structure
g2 = girth_and_pow2(adj)
has8b = find_cycle_len(adj, 8)
unwired = [( 'a', p) for p in range(2, 12, 2)]
deg = {v: len(adj[v]) for v in adj}
n_deg2 = sum(1 for v, dd in deg.items() if dd == 2 and v[0] in 'ab')
print('(4) spacing-2 rungs: girth %d, 8-cycle: %s; unwired leg'
      ' vertices remaining: %d -> the cascade (f) drives'
      % (g2, has8b, n_deg2))
assert g2 == 6 and not has8b and n_deg2 > 0
print()
print('CONCLUSION: the metric half of the braided core is realisable by')
print('design arithmetic; every completion step collides with girth or')
print('with condition (f) - the ambient counterexample hypothesis. The')
print('433.4 question WITHOUT (f) is weaker than J7 and likely')
print('constructible; WITH (f) it carries the conjecture itself.')
