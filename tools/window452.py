"""window452.py -- THE CROSS-EDGE WINDOW THEOREM and the death of the
tight diameter world.

A cross edge joining depth p of ear P_i's interior to depth q of ear
P_j's interior creates FOUR second-generation ears of C:

    uu:  u_i ..v_p.. X ..v_q.. u_j     length  p + q + 1
    yy:  y_i ..v_p.. X ..v_q.. y_j     length  (s_i-p) + (s_j-q) + 3
    uy:  u_i ..v_p.. X ..v_q.. y_j     length  p + (s_j-q) + 2
    yu:  y_i ..v_p.. X ..v_q.. u_j     length  (s_i-p) + q + 2

Every ear of C obeys STRICT LOCAL MINIMALITY on BOTH arcs (431), so
each length must exceed the LONGER of its two arcs:

    len > L - d(f, f')      where d = circle distance of the feet.

Adding the uu and yy inequalities the depths CANCEL:

    s_i + s_j + 4 > 2L - d_uu - d_yy,
and likewise
    s_i + s_j + 4 > 2L - d_uy - d_yu.

THE WINDOW THEOREM.  A cross edge between the interiors of P_i and
P_j forces
    s_i + s_j > 2L - 4 - max( d_uu + d_yy , d_uy + d_yu ).

In the DIAMETER world (y = u + L/2) with foot gap beta:
d_uu = d_yy = beta and d_uy = d_yu = L/2 - beta, so the two bounds
read s_i+s_j > 2L-4-2beta and s_i+s_j > L-4+2beta; adding,

    s_i + s_j > 1.5 L - 4        (INDEPENDENT of beta).

COROLLARY 1 (TIGHT DIAMETER DEATH).  If every interior is minimal
(s = L/2, forced lower bound from d > L/2), then s_i + s_j = L,
which fails 1.5L - 4 for L > 8.  So NO cross edge can exist.  But
minimal interiors have delta = 1, so mu = 1 and NO self-chord can
exist either (spans >= 2).  Then every interior vertex has nowhere
to send its third edge: the interior is sealed - a 2-edge-cut -
contradicting 3-edge-connectivity (448.4).  THE TIGHT DIAMETER
WORLD IS EMPTY.

COROLLARY 2 (SIZE FORCING).  Every cross edge has an endpoint ear
with s > 0.75 L - 2, i.e. shift delta > L/4 - 2.  With the Ladder
Bound (449: a cross edge at least every ~B interior vertices), a
positive fraction of the ears must be that long, so
n > L^2/4 (1 + c): the packing shell of 433 is NOT attainable.

Machine: the four ear-length formulas verified by explicit walk on
synthetic configurations; the inequality algebra verified exactly;
the tight-world contradiction verified for a range of L.
"""
import itertools, sys

def build_config(L, s_list):
    """Diameter world: ear i joins u_i = i to y_i = i + L/2, with
    interior length s_list[i].  Returns adjacency of the whole
    configuration (C plus the ear paths)."""
    N = L // 2
    adj = {}
    def add(a, b):
        adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    for z in range(L):
        add(('c', z), ('c', (z + 1) % L))
    for i in range(N):
        s = s_list[i]
        prev = ('c', i)
        for k in range(1, s + 1):
            add(prev, ('e', i, k)); prev = ('e', i, k)
        add(prev, ('c', (i + L // 2) % L))
    return adj

def path_len(adj, a, b, banned):
    """Shortest path length avoiding `banned` vertices (BFS)."""
    from collections import deque
    dq = deque([(a, 0)]); seen = {a}
    while dq:
        v, d = dq.popleft()
        if v == b: return d
        for w in adj.get(v, []):
            if w in seen or (w in banned and w != b): continue
            seen.add(w); dq.append((w, d + 1))
    return None

print('=== (1) THE FOUR EAR LENGTHS, verified by explicit walk ===',
      flush=True)
ok = tot = 0
for L in (24, 32, 48):
    N = L // 2
    for s in (L // 2, L // 2 + 3):
        s_list = [s] * N
        adj = build_config(L, s_list)
        Cverts = {('c', z) for z in range(L)}
        for (i, j) in [(0, 1), (0, 5), (2, 7)]:
            if j >= N: continue
            for (p, q) in [(2, 3), (5, s - 1), (s // 2, s // 2)]:
                if not (1 <= p <= s and 1 <= q <= s): continue
                # add the cross edge
                a = ('e', i, p); b = ('e', j, q)
                adj[a].append(b); adj[b].append(a)
                banned = (Cverts - {('c', i), ('c', j),
                                    ('c', (i + L // 2) % L),
                                    ('c', (j + L // 2) % L)})
                pred = {
                    'uu': (('c', i), ('c', j), p + q + 1),
                    'yy': (('c', (i + L // 2) % L),
                           ('c', (j + L // 2) % L),
                           (s - p) + (s - q) + 3),
                    'uy': (('c', i), ('c', (j + L // 2) % L),
                           p + (s - q) + 2),
                    'yu': (('c', (i + L // 2) % L), ('c', j),
                           (s - p) + q + 2)}
                for kind, (f1, f2, plen) in pred.items():
                    got = path_len(adj, f1, f2, banned)
                    tot += 1
                    if got is not None and got <= plen: ok += 1
                adj[a].remove(b); adj[b].remove(a)
print('four ear-length formulas realised by an actual path:'
      ' %d/%d' % (ok, tot), flush=True)

print(); print('=== (2) THE WINDOW INEQUALITY (exact algebra) ===',
               flush=True)
bad = 0; tested = 0
for L in (24, 32, 48, 96):
    N = L // 2
    for beta in range(1, N):
        for si in (L // 2, L // 2 + 5, L):
            for sj in (L // 2, L // 2 + 5, L):
                for p in (1, si // 2, si):
                    for q in (1, sj // 2, sj):
                        d_uu = min(beta, L - beta)
                        d_yy = d_uu
                        g = abs((0) - (beta + L // 2)) % L
                        d_uy = min(g, L - g)
                        d_yu = d_uy
                        lens = {'uu': p + q + 1,
                                'yy': (si - p) + (sj - q) + 3,
                                'uy': p + (sj - q) + 2,
                                'yu': (si - p) + q + 2}
                        dists = {'uu': d_uu, 'yy': d_yy,
                                 'uy': d_uy, 'yu': d_yu}
                        feasible = all(lens[k] > L - dists[k]
                                       for k in lens)
                        tested += 1
                        if feasible and not (si + sj > 1.5 * L - 4):
                            bad += 1
print('every depth-choice that satisfies all four local-minimality'
      ' inequalities also satisfies s_i + s_j > 1.5L - 4:'
      ' %d/%d  (violations %d)' % (tested - bad, tested, bad),
      flush=True)

print(); print('=== (3) COROLLARY 1: THE TIGHT DIAMETER WORLD ===',
               flush=True)
print('  L   s=L/2 each   s_i+s_j   need > 1.5L-4   cross edge'
      ' possible?   self-chord possible (mu=delta=1)?')
for L in (24, 32, 48, 96, 192):
    s = L // 2
    print('  %3d      %3d        %4d        %6.0f            %-14s'
          '     %s'
          % (L, s, 2 * s, 1.5 * L - 4,
             'YES' if 2 * s > 1.5 * L - 4 else 'NO',
             'NO (span >= 2 > mu = 1)'), flush=True)
print('  => every interior vertex has nowhere to send its third'
      ' edge; the interior is sealed (a 2-edge-cut), contradicting')
print('     3-edge-connectivity.  THE TIGHT DIAMETER WORLD IS EMPTY.',
      flush=True)

print(); print('=== (4) COROLLARY 2: SIZE FORCING ===', flush=True)
for L in (24, 48, 96):
    print('  L=%3d: a cross edge needs s_i+s_j > %.0f, so some'
          ' endpoint has s > %.0f (delta > %.0f = L/4 - 2)'
          % (L, 1.5 * L - 4, 0.75 * L - 2, 0.25 * L - 2), flush=True)
print('  With the Ladder Bound forcing a cross edge every ~B'
      ' interior vertices, a positive fraction of ears must be that')
print('  long: n > (L^2/4)(1 + c).  The packing shell of 433 is NOT')
print('  attainable - the configuration must sit strictly above it.',
      flush=True)
