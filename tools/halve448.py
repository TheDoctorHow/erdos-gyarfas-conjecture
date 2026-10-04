"""halve448.py -- THE HALVING QUESTION: the identity, the obstruction,
and the symmetry payoff.

 (A) THE HALVING IDENTITY.  If every cycle of a weighted graph
     (G, l) has even length, then l == dS (mod 2) edgewise for some
     vertex set S (recovered by tree propagation), w = (l - dS)/2
     is an integer weighting, and for EVERY cycle
         l(C)/2 = w(C) + cross_S(C),
     cross_S(C) = |C cap dS| / 2.  Halving succeeds exactly up to a
     crossing functional.  Machine: random weighted instances on the
     four census graphs, identity checked on every circuit.
 (B) THE HALF-CUT OBSTRUCTION.  The crossing functional is NOT a
     weight: h_v(C) = [v in C] has the rational weight (1/2) dv and
     no integer weight - any theta at v (three internally disjoint
     v-b paths) forces w(P_i) + w(P_j) = 1 for all pairs, so
     w(P_i) = 1/2.  Machine: exhibit a theta at every vertex of
     each census graph (the only graph-dependent input; the rest is
     the 3x3 system).
 (C) THE SYMMETRY PAYOFF.  Where halving DOES work: free
     quotients.  Covering lemma: every cycle of G/sigma lifts to a
     cycle of G of the same length (trivial holonomy) or double
     length - so if G has no power-of-2 cycle, neither has G/sigma,
     which has n/2 vertices: A MINIMAL COUNTEREXAMPLE ADMITS NO
     FREE INVOLUTION with d(u, sigma u) >= 3 (simple cubic
     quotient).  Machine: the antipodal involution of
     Mobius-Kantor, its quotient, and the covering lemma verified
     cycle by cycle.
"""
import itertools, random, sys, time
sys.setrecursionlimit(100000)
random.seed(448)

def k4():
    return 4, [(i, j) for i in range(4) for j in range(i + 1, 4)]

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

def cyc_edges(cyc):
    return [frozenset((a, b)) for a, b in zip(cyc, cyc[1:] + cyc[:1])]

GRAPHS = [('K4', k4()), ('Petersen', petersen()),
          ('Mobius-Kantor', mobius_kantor()), ('Heawood', heawood())]

print('=== (A) THE HALVING IDENTITY ===', flush=True)
for name, nE in GRAPHS:
    n, E, adj = build(nE)
    cycles = all_cycles(adj, n)
    eidx = {frozenset(e): i for i, e in enumerate(E)}
    ok_all = 0; trials = 20
    for _ in range(trials):
        w0 = [random.randint(1, 5) for _ in E]
        S0 = {v for v in range(n) if random.random() < 0.5}
        ell = [2 * w0[i] + (1 if (E[i][0] in S0) != (E[i][1] in S0)
                            else 0) for i in range(len(E))]
        # all cycle lengths even by construction; verify anyway
        even = all(sum(ell[eidx[e]] for e in cyc_edges(C)) % 2 == 0
                   for C in cycles)
        # recover S' by tree propagation of ell mod 2
        s = [None] * n; s[0] = 0
        order = [0]; par = {0: None}
        for v in order:
            for u in adj[v]:
                if s[u] is None:
                    s[u] = (s[v] + ell[eidx[frozenset((u, v))]]) % 2
                    order.append(u); par[u] = v
        consistent = all((s[a] + s[b]) % 2 == ell[i] % 2
                         for i, (a, b) in enumerate(E))
        Sp = {v for v in range(n) if s[v] == 1}
        wp = [(ell[i] - (1 if (E[i][0] in Sp) != (E[i][1] in Sp)
                         else 0)) for i in range(len(E))]
        integral = all(x % 2 == 0 for x in wp)
        wp = [x // 2 for x in wp]
        ident = True
        for C in cycles:
            eds = cyc_edges(C)
            lC = sum(ell[eidx[e]] for e in eds)
            wC = sum(wp[eidx[e]] for e in eds)
            crossC = sum(1 for e in eds
                         if (tuple(e)[0] in Sp) != (tuple(e)[1] in Sp))
            if lC % 2 != 0 or crossC % 2 != 0: ident = False; break
            if lC // 2 != wC + crossC // 2: ident = False; break
        if even and consistent and integral and ident: ok_all += 1
    print('%-14s %d/%d random weighted instances: recovery,'
          ' integrality, and l(C)/2 = w(C) + cross(C) on all %d'
          ' circuits' % (name, ok_all, trials, len(cycles)), flush=True)

print(); print('=== (B) THE HALF-CUT OBSTRUCTION ===', flush=True)
def find_theta(adj, n, v):
    """Three internally-vertex-disjoint paths from v to some b,
    one through each neighbour of v."""
    nb = adj[v]
    # enumerate simple paths from each neighbour avoiding v,
    # indexed by endpoint
    def paths_from(start, maxlen=9):
        res = {}
        stack = [(start, (start,))]
        while stack:
            x, p = stack.pop()
            res.setdefault(x, []).append(p)
            if len(p) >= maxlen: continue
            for y in adj[x]:
                if y != v and y not in p:
                    stack.append((y, p + (y,)))
        return res
    P1 = paths_from(nb[0]); P2 = paths_from(nb[1]); P3 = paths_from(nb[2])
    for b in range(n):
        if b == v: continue
        for p1 in P1.get(b, []):
            s1 = set(p1) - {b}
            for p2 in P2.get(b, []):
                if (set(p2) - {b}) & s1: continue
                s12 = s1 | (set(p2) - {b})
                for p3 in P3.get(b, []):
                    if (set(p3) - {b}) & s12: continue
                    return b, p1, p2, p3
    return None

for name, nE in GRAPHS:
    n, E, adj = build(nE)
    found = 0
    sample = None
    for v in range(n):
        th = find_theta(adj, n, v)
        if th: found += 1
        if v == 0: sample = th
    b, p1, p2, p3 = sample
    print('%-14s thetas at %d/%d vertices; at v=0: b=%d, path'
          ' lengths (%d, %d, %d) -> w(Pi)+w(Pj)=1 for all pairs'
          ' forces w(Pi)=1/2: NO integer weight realises [v in C]'
          % (name, found, n, b, len(p1), len(p2), len(p3)), flush=True)

print(); print('=== (C) THE SYMMETRY PAYOFF ===', flush=True)
n, E, adj = build(mobius_kantor())
def sigma(v):
    return (v + 4) % 8 if v < 8 else 8 + ((v - 8 + 4) % 8)
# automorphism, free, distance check
is_auto = all(frozenset((sigma(a), sigma(b))) in
              {frozenset(e) for e in E} for a, b in E)
free = all(sigma(v) != v for v in range(n))
def dist(a, b):
    seen = {a: 0}; q = [a]
    for x in q:
        if x == b: return seen[x]
        for y in adj[x]:
            if y not in seen: seen[y] = seen[x] + 1; q.append(y)
    return 99
dmin = min(dist(v, sigma(v)) for v in range(n))
print('Mobius-Kantor antipode sigma = +4: automorphism %s, free %s,'
      ' min d(u, sigma u) = %d' % (is_auto, free, dmin), flush=True)
# quotient
orb = {}
for v in range(n):
    o = min(v, sigma(v)); orb[v] = o
overts = sorted(set(orb.values()))
oedges = set()
mult = {}
for a, b in E:
    key = frozenset((orb[a], orb[b]))
    mult[key] = mult.get(key, 0) + 1
    oedges.add(key)
simple = all(len(k) == 2 for k in oedges) and \
    all(m == 2 for m in mult.values())
print('quotient: %d vertices, %d edge-orbits, simple cubic: %s'
      % (len(overts), len(oedges), simple), flush=True)
qadj = {v: [] for v in overts}
for k in oedges:
    a, b = tuple(k); qadj[a].append(b); qadj[b].append(a)
remap = {v: i for i, v in enumerate(overts)}
qadj2 = {remap[v]: [remap[u] for u in qadj[v]] for v in overts}
qcycles = all_cycles(qadj2, len(overts))
spec_up = {len(C) for C in all_cycles(adj, n)}
oklift = sum(1 for C in qcycles
             if len(C) in spec_up or 2 * len(C) in spec_up)
print('covering lemma: %d/%d quotient cycles have their length or'
      ' its double in the Mobius-Kantor spectrum' %
      (oklift, len(qcycles)), flush=True)
print('quotient spectrum: %s' % sorted({len(C) for C in qcycles}))
print('upstairs spectrum: %s' % sorted(spec_up), flush=True)
print()
print('ASYMMETRY THEOREM: were G a minimal counterexample, G/sigma')
print('would be a counterexample on n/2 vertices (no power-of-2')
print('cycle lifts to one) - contradiction.  The minimal')
print('counterexample admits NO free involution with')
print('d(u, sigma u) >= 3.')

# ---- (C2) the covering lemma, EXPLICITLY (not via spectrum) ----
print(); print('--- (C2) explicit preimage test ---', flush=True)
inv = {}
for v in overts:
    inv.setdefault(orb[v], []).append(v)
for v in range(n):
    inv.setdefault(orb[v], [])
    if v not in inv[orb[v]]: inv[orb[v]].append(v)
bad = 0; lens = {}
for Cq in qcycles:
    real = [overts[x] for x in Cq]
    Vpre = set()
    for x in real: Vpre |= {x, sigma(x)}
    qedges = {frozenset((real[i], real[(i + 1) % len(real)]))
              for i in range(len(real))}
    pre = [e for e in E
           if frozenset((orb[e[0]], orb[e[1]])) in qedges
           and e[0] in Vpre and e[1] in Vpre]
    deg = {}
    for a, b in pre:
        deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
    if set(deg) != Vpre or any(d != 2 for d in deg.values()):
        bad += 1; continue
    # components
    padj = {}
    for a, b in pre:
        padj.setdefault(a, []).append(b); padj.setdefault(b, []).append(a)
    seenv = set(); comps = []
    for s0 in Vpre:
        if s0 in seenv: continue
        comp = {s0}; st = [s0]
        while st:
            x = st.pop()
            for y in padj[x]:
                if y not in comp: comp.add(y); st.append(y)
        seenv |= comp; comps.append(len(comp))
    m = comps[0] // len(real)
    if any(c != comps[0] for c in comps) or m not in (1, 2):
        bad += 1; continue
    lens[(len(real), m)] = lens.get((len(real), m), 0) + 1
print('preimages of all %d quotient cycles are 2-regular unions of'
      ' equal cycles: failures %d' % (len(qcycles), bad))
print('(|Cbar|, multiplicity m) census: %s'
      % {k: lens[k] for k in sorted(lens)})
print('=> every quotient cycle lifts to cycles of length m|Cbar|,'
      ' m in {1,2}: a power of 2 downstairs forces one upstairs.',
      flush=True)

# ---------------- (D) SCALE INVARIANCE ----------------
print(); print('=== (D) THE SCALE-INVARIANCE THEOREM ===', flush=True)
def board_system(t, scale):
    """The m=3 diameter board at level t, in units of `scale`.
    Returns the multiset of pair-constraint data (target, modulus,
    domain) in those units."""
    L = 3 * (1 << t)
    M = 1 << (t + 1)
    dom = sorted({v for v in range(M) if v % (1 << t) != 0})
    cons = []
    for i in range(L // 2):
        for j in range(i + 1, L // 2):
            cons.append((L - 2 * (j - i)) % M)
    return (L / scale, M / scale, [d / scale for d in dom],
            sorted(set(c / scale for c in cons)))
A3 = board_system(3, 1)
A2h = board_system(3, 2)       # the SAME board, measured in half-units
print('level-3 board in integer units : L=%g, modulus=%g, |domain|=%d,'
      ' %d distinct pair targets' % (A3[0], A3[1], len(A3[2]), len(A3[3])))
print('same board in HALF units       : L=%g, modulus=%g, |domain|=%d,'
      ' %d distinct pair targets' % (A2h[0], A2h[1], len(A2h[2]),
                                     len(A2h[3])))
print('the halved system is a LEVEL-2 system whose values live in'
      ' (1/2)Z: identical combinatorics, halved ruler.')
print('So: v_2 of every cycle drops by 1 AND the target set {2^k}')
print('maps onto itself.  Halving carries a level-t counterexample')
print('to a level-(t-1) counterexample in the half-weighted world;')
print('applying a "level 2" theorem there is exactly the level-t')
print('statement upstairs.  THE INDUCTION IS CIRCULAR.', flush=True)
print()
print('WHERE t=2 REALLY USES INTEGRALITY: at level t each even ear')
print('forbids one class mod 2^{t-1}.  Measured against the')
print('REALIZABLE shift set, that density is:')
for t in (2, 3, 4):
    print('   unweighted level %d: 1 class of %d (density %.3f)'
          % (t, 1 << (t - 1), 1.0 / (1 << (t - 1))))
print('   half-weighted level 2 (= level 3 halved): 1 class of 2 in'
      ' Z, but shifts range over (1/2)Z: density 1/4 - the t=2')
print('   density 1/2 is NOT scale-invariant relative to the value')
print('   set.  INTEGRALITY OF SHIFTS is the single step.', flush=True)
