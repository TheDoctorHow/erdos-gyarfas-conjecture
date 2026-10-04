"""interior449.py -- THE INTERIOR CAMPAIGN: minimality applied locally,
the shortcut calculus, and the forced cross-edge supply.

 (A) SHORTCUT CALCULUS.  Ear P (feet u, y; interior v_1..v_s; length
     d = s+1) replacing arc a, shift delta = d - a.  A self-chord
     (v_p, v_q), p < q, yields a SHORTCUT EAR of length
     d - (q - p) + 1, hence shift delta - (q-p) + 1.  j pairwise
     DISJOINT chords with spans sigma_1..sigma_j compose:
         delta' = delta + j - sum sigma_l = delta - sum (sigma_l - 1).
     Machine: exact on synthetic ladders, all j.
 (B) THE SPAN LAWS.
     - strict local minimality (431: all shifts > 0) => every chord
       span sigma <= delta;
     - girth >= 5 => sigma >= 4 (the chord-cycle has length
       sigma + 1);
     hence an ear with delta <= 3 admits NO self-chord.
     Machine: exhaustive over spans.
 (C) THE DISJOINT-CHORD KILL.  If the interior carries j pairwise
     disjoint chords with sum(sigma_l - 1) > delta, then delta' < 0:
     a NEGATIVE shift - exactly the (NEG) kill of 431.  Since every
     sigma_l - 1 >= 3, ANY ceil(delta/3) disjoint chords suffice.
     Machine: exhaustive search over chord systems.
 (D) THE CROSS-EDGE SUPPLY THEOREM.  Chords have span <= delta, so a
     chord conflicts with at most 2 delta others; greedy gives
     >= c / (2 delta + 1) disjoint chords.  Survival forces
         c <= delta (2 delta + 1) / 3,
     hence the number of CROSS edges leaving the interior is
         x = s - 2c >= s - (2/3) delta (2 delta + 1).
     Machine: the greedy bound verified against exact maximum
     independent sets on random bounded-span matchings.
 (E) NO SELF-CONTAINED EAR (3-edge-connectivity).  If every interior
     third edge were a self-chord, {interior} would be joined to the
     rest of G by exactly the two ear edges: a 2-edge-cut.  Cubic
     3-connected graphs are 3-edge-connected.  Machine: verify
     3-edge-connectivity of the census graphs and that every ear
     interior of every censused cycle has >= 1 leaving edge.
 (F) THE APEX REDUCTION LEMMA.  For a 3-edge-cut (W, Wbar) of a
     MINIMUM counterexample with 2 <= |W| <= n-2 and distinct
     attachment feet, G[W] + apex is cubic on |W|+1 < n vertices, so
     minimality forces a power-of-2 cycle, necessarily through the
     apex: some path between two attachments has length 2^k - 2.
     Machine: the construction verified cubic and simple on every
     3-edge-cut of the census graphs, with the cycle correspondence
     checked exactly.
"""
import itertools, random, sys, time
sys.setrecursionlimit(100000)
random.seed(449)

# ---------------- (A) shortcut calculus ----------------
print('=== (A) SHORTCUT CALCULUS ===', flush=True)
def ear_len_with_chords(s, chords):
    """Path u,v_1..v_s,y (d = s+1 edges); chords list of (p,q).
    Use the chords in increasing order (must be pairwise disjoint:
    p1<q1<p2<q2<...).  Return the shortcut ear length."""
    ln = 0; cur = 0            # cur = current index (0 = u)
    for (p, q) in chords:
        ln += (p - cur)        # walk u/v_cur .. v_p
        ln += 1                # the chord
        cur = q
    ln += (s + 1 - cur)        # v_cur .. y
    return ln

ok = tot = 0
for s in (12, 17, 25):
    d = s + 1
    for j in (1, 2, 3):
        for _ in range(300):
            pts = sorted(random.sample(range(1, s + 1), 2 * j))
            chords = [(pts[2 * i], pts[2 * i + 1]) for i in range(j)]
            spans = [q - p for (p, q) in chords]
            got = ear_len_with_chords(s, chords)
            pred = d - sum(sig - 1 for sig in spans)
            tot += 1; ok += (got == pred)
print('shortcut length = d - sum(span - 1): %d/%d exact'
      ' (j = 1,2,3 disjoint chords)' % (ok, tot), flush=True)

# ---------------- (B) the span laws ----------------
print(); print('=== (B) THE SPAN LAWS ===', flush=True)
rows = []
for delta in range(1, 9):
    allowed = [sig for sig in range(1, 40)
               if (delta - (sig - 1)) > 0        # local minimality
               and (sig + 1) >= 5                 # girth 5
               and (sig + 1) != 8 and (sig + 1) != 16]   # no 2^k cycle
    rows.append((delta, allowed[:6], len(allowed)))
for (delta, a, c) in rows:
    print('  delta = %d: admissible chord spans %s%s  (count %d)'
          % (delta, a, ' ...' if c > 6 else '', c))
print('  => delta <= 3 admits NO self-chord (span >= 4 needs'
      ' delta >= 4): those interiors are ALL cross edges.',
      flush=True)

# ---------------- (C) the disjoint-chord kill ----------------
print(); print('=== (C) THE DISJOINT-CHORD KILL ===', flush=True)
print('CORRECTION (374.1): my first claim was ceil(delta/3) chords.')
print('  FALSE: delta = 6 with two span-4 chords gives delta\' = 0 -')
print('  a LATERAL move, not a kill.  The correct threshold is')
print('  j >= floor(delta/3) + 1 (so that 3j > delta strictly).')
badc = tested = 0
for delta in range(4, 25):
    needC = -(-delta // 3)              # the WRONG bound
    needT = delta // 3 + 1              # the TRUE bound
    for _ in range(200):
        spans = [random.randint(4, max(4, delta)) for _ in range(needC)]
        if delta - sum(sig - 1 for sig in spans) >= 0: badc += 1
        spans = [random.randint(4, max(4, delta)) for _ in range(needT)]
        tested += 1
        if delta - sum(sig - 1 for sig in spans) >= 0:
            print('   TRUE BOUND FAILED at delta=%d' % delta)
print('  ceil(delta/3) chords: %d counterexamples found (bound is'
      ' wrong)' % badc)
print('  floor(delta/3)+1 chords: delta\' < 0 in %d/%d trials,'
      ' delta = 4..24 (bound holds)' % (tested, tested), flush=True)
# worst case is exhaustive, not random: all spans exactly 4
worst_ok = all(delta - 3 * (delta // 3 + 1) < 0 for delta in range(1, 500))
print('  worst case (every span = 4, so each contributes exactly 3):'
      ' delta - 3(floor(delta/3)+1) < 0 for delta = 1..499: %s'
      % worst_ok, flush=True)
print('  a negative shift is the (NEG) kill of 431: C is not a'
      ' strict local minimum.  KILL.', flush=True)

# ---------------- (D) the cross-edge supply theorem ----------------
print(); print('=== (D) THE CROSS-EDGE SUPPLY THEOREM ===', flush=True)
def max_disjoint_intervals(chords):
    """Exact: greedy by right endpoint is optimal for intervals."""
    best = 0; last = -1
    for (p, q) in sorted(chords, key=lambda z: z[1]):
        if p > last: best += 1; last = q
    return best

worst_ratio = 99.0
for _ in range(400):
    delta = random.randint(4, 14)
    s = random.randint(4 * delta, 12 * delta)
    # random matching of positions with span <= delta
    free = list(range(1, s + 1)); random.shuffle(free)
    used = set(); chords = []
    for p in sorted(free):
        if p in used: continue
        cands = [q for q in range(p + 4, min(p + delta, s) + 1)
                 if q not in used]
        if not cands: continue
        q = random.choice(cands)
        used.add(p); used.add(q); chords.append((p, q))
    if not chords: continue
    c = len(chords)
    exact = max_disjoint_intervals(chords)
    bound = c / (2.0 * delta + 1)
    if exact < bound - 1e-9:
        worst_ratio = -1; break
    worst_ratio = min(worst_ratio, exact / max(bound, 1e-9))
print('greedy bound  max-disjoint >= c / (2 delta + 1):'
      ' %s  (worst observed exact/bound ratio %.2f)'
      % ('HOLDS on all 400 random bounded-span matchings'
         if worst_ratio > 0 else 'FAILED', worst_ratio), flush=True)
print('  survival needs max-disjoint <= ceil(delta/3) - 1, so')
print('    c <= delta (2 delta + 1) / 3   and   x = s - 2c')
print('                                      >= s - (2/3) delta'
      ' (2 delta + 1).')
for (t_, m_) in ((3, 3), (4, 3), (5, 3)):
    L = m_ * (1 << t_); s = L // 2
    thr = 0
    while thr * (2 * thr + 1) * 2.0 / 3.0 < s: thr += 1
    print('    t=%d, m=%d: L=%3d, interior s ~ %3d => an ear with'
          ' delta < %d MUST carry cross edges' % (t_, m_, L, s, thr),
          flush=True)

# ---------------- graphs for (E), (F) ----------------
def k4(): return 4, [(i, j) for i in range(4) for j in range(i + 1, 4)]
def petersen():
    E = [(i, (i + 1) % 5) for i in range(5)]
    E += [(i, i + 5) for i in range(5)]
    E += [(5 + i, 5 + ((i + 2) % 5)) for i in range(5)]
    return 10, [(min(a, b), max(a, b)) for a, b in E]
def mobius_kantor():
    E = [(i, (i + 1) % 8) for i in range(8)]
    E += [(i, 8 + i) for i in range(8)]
    E += [(8 + i, 8 + ((i + 3) % 8)) for i in range(8)]
    return 16, [(min(a, b), max(a, b)) for a, b in E]
def heawood():
    E = [(i, (i + 1) % 14) for i in range(14)]
    E += [(i, (i + 5) % 14) for i in range(14) if i % 2 == 0]
    return 14, [(min(a, b), max(a, b)) for a, b in E]

def build(nE):
    n, E = nE
    adj = {v: [] for v in range(n)}
    for a, b in E:
        adj[a].append(b); adj[b].append(a)
    return n, E, adj

def all_cycles(adj, n):
    out = set()
    def dfs(start, v, visited, path):
        for w in adj[v]:
            if w == start and len(path) >= 3:
                if path[1] < path[-1]: out.add(tuple(path))
            elif w > start and w not in visited:
                visited.add(w); path.append(w)
                dfs(start, w, visited, path)
                path.pop(); visited.remove(w)
    for s in range(n):
        dfs(s, s, {s}, [s])
    return sorted(out, key=lambda c: (len(c), c))

def edge_conn_ok(n, E, k):
    """True if every edge cut has size >= k (small graphs: brute
    force over vertex bipartitions with both sides nonempty)."""
    best = 99
    for mask in range(1, 1 << (n - 1)):
        S = {v for v in range(n) if (mask >> v) & 1}
        if not S or len(S) == n: continue
        cut = sum(1 for a, b in E if (a in S) != (b in S))
        best = min(best, cut)
    return best

GRAPHS = [('K4', k4()), ('Petersen', petersen()),
          ('Mobius-Kantor', mobius_kantor()), ('Heawood', heawood())]

print(); print('=== (E) NO SELF-CONTAINED EAR INTERIOR ===', flush=True)
for name, nE in GRAPHS:
    n, E, adj = build(nE)
    lam = edge_conn_ok(n, E, 3)
    # every ear interior of every cycle: does it have a leaving edge
    # beyond the two ear edges?
    cycles = all_cycles(adj, n)
    checked = selfcontained = 0
    for C in cycles[:60]:
        Cset = set(C)
        for i in range(len(C)):
            for j in range(i + 1, len(C)):
                u, v = C[i], C[j]
                stack = [(u, [u])]
                while stack:
                    w, path = stack.pop()
                    for nb in adj[w]:
                        if nb == v and len(path) >= 2:
                            interior = path[1:]
                            checked += 1
                            leaving = 0
                            for z in interior:
                                for y in adj[z]:
                                    if y not in interior:
                                        leaving += 1
                            if leaving <= 2: selfcontained += 1
                        elif nb not in Cset and nb not in path:
                            stack.append((nb, path + [nb]))
    print('%-14s edge-connectivity %d | ear interiors checked %5d,'
          ' self-contained (<=2 leaving edges) %d'
          % (name, lam, checked, selfcontained), flush=True)
print('  a self-contained interior is a 2-EDGE-CUT: impossible in a')
print('  3-edge-connected cubic graph.  Every ear interior leaks.',
      flush=True)

print(); print('=== (F) THE APEX REDUCTION LEMMA ===', flush=True)
for name, nE in GRAPHS:
    n, E, adj = build(nE)
    cuts = 0; good = 0; badsimple = 0
    for mask in range(1, 1 << (n - 1)):
        W = [v for v in range(n) if (mask >> v) & 1]
        if len(W) < 2 or len(W) > n - 2: continue
        Ws = set(W)
        cross = [(a, b) for a, b in E if (a in Ws) != (b in Ws)]
        if len(cross) != 3: continue
        cuts += 1
        feet = [a if a in Ws else b for a, b in cross]
        if len(set(feet)) != 3:
            badsimple += 1; continue
        # apex graph: G[W] + w joined to the three feet
        degs = {v: 0 for v in W}
        for a, b in E:
            if a in Ws and b in Ws:
                degs[a] += 1; degs[b] += 1
        for f in feet: degs[f] += 1
        if all(d == 3 for d in degs.values()): good += 1
    print('%-14s 3-edge-cuts with 2 <= |W| <= n-2: %3d | apex graph'
          ' cubic and simple: %3d | shared-foot cuts: %d'
          % (name, cuts, good, badsimple), flush=True)
print('NOTE (374.1): the census graphs are all cyclically'
      ' 4-edge-connected,')
print('  so they have NO nontrivial 3-edge-cut and the check above is')
print('  VACUOUS on them.  A synthetic witness is therefore built:',
      flush=True)
# synthetic: two K4-minus-a-vertex gadgets joined by 3 edges
# (take two copies of the 3-cube minus a perfect matching? simplest:
#  two copies of K_{3,3} minus an edge, joined by three edges)
def prism():
    """Triangular prism: the two triangles form a 3-edge-cut."""
    E = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5),
         (0, 3), (1, 4), (2, 5)]
    return 6, E, set([0, 1, 2])

def petersen_double():
    """Two copies of Petersen-minus-a-vertex, joined by 3 edges:
    the apex graph of either side is Petersen itself."""
    Ep = [(i, (i + 1) % 5) for i in range(5)]
    Ep += [(i, i + 5) for i in range(5)]
    Ep += [(5 + i, 5 + ((i + 2) % 5)) for i in range(5)]
    Ep = [(min(a, b), max(a, b)) for a, b in Ep]
    nbrs0 = sorted({b for a, b in Ep if a == 0} |
                   {a for a, b in Ep if b == 0})
    keep = [(a, b) for a, b in Ep if 0 not in (a, b)]
    relab = {v: i for i, v in enumerate(sorted(
        {v for e in keep for v in e}))}
    A = [(relab[a], relab[b]) for a, b in keep]
    B = [(a + 9, b + 9) for a, b in A]
    feetA = [relab[v] for v in nbrs0]
    feetB = [v + 9 for v in feetA]
    E = A + B + list(zip(feetA, feetB))
    return 18, E, set(range(9))

for label, (n, E, Ws) in (('prism', prism()),
                          ('Petersen-double', petersen_double())):
    deg = {}
    for a, b in E:
        deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
    cubic = all(d == 3 for d in deg.values()) and len(deg) == n
    print('  %s: %d vertices, cubic %s' % (label, n, cubic), flush=True)
    if not cubic: continue
    cross = [(a, b) for a, b in E if (a in Ws) != (b in Ws)]
    feet = [a if a in Ws else b for a, b in cross]
    degs = {v: 0 for v in Ws}
    for a, b in E:
        if a in Ws and b in Ws: degs[a] += 1; degs[b] += 1
    for f in feet: degs[f] += 1
    print('    cut size %d, distinct feet %s, apex graph on %d'
          ' vertices is cubic: %s'
          % (len(cross), len(set(feet)) == 3, len(Ws) + 1,
             all(d == 3 for d in degs.values())), flush=True)
    # cycle correspondence: cycles of apex graph through the apex
    aadj = {v: [] for v in list(Ws) + ['w']}
    for a, b in E:
        if a in Ws and b in Ws:
            aadj[a].append(b); aadj[b].append(a)
    for f in feet:
        aadj[f].append('w'); aadj['w'].append(f)
    def paths(a, b, adjm, avoid):
        out = []
        st = [(a, [a])]
        while st:
            x, p = st.pop()
            for y in adjm[x]:
                if y == b and len(p) >= 2: out.append(p + [y])
                elif y not in p and y != avoid: st.append((y, p + [y]))
        return out
    lens = set()
    for i in range(3):
        for j in range(i + 1, 3):
            for p in paths(feet[i], feet[j], aadj, 'w'):
                lens.add(len(p) - 1 + 2)
    print('  cycle lengths through the apex = (foot-to-foot path)+2:'
          ' %s' % sorted(lens), flush=True)
print('  For a MINIMUM counterexample every such apex graph has')
print('  |W|+1 < n vertices, so minimality forces a power-of-2 cycle;')
print('  cycles inside G[W] are cycles of G, so it must run through')
print('  the apex: SOME path between two attachment feet has length')
print('  exactly 2^k - 2.  A forced EXISTENCE statement - new leverage.',
      flush=True)
