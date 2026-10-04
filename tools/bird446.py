"""bird446.py -- THE BIRD PROGRAMME I: the spectrum reformulation,
the One-Step Theorem, and the island measurements.

 (A) TWO CUBIC FACTS.  In a cubic graph: (i) two circuits sharing a
     vertex share an edge at that vertex (two 2-subsets of a 3-set
     meet); (ii) at any vertex of C, a circuit C' through it uses at
     least one C-edge.  Machine: all circuit pairs, four graphs.
 (B) ONE-STEP THEOREM.  For any circuits C, C' of a cubic graph with
     V(C) cap V(C') nonempty: C' decomposes as (kept arcs of C) +
     (ears P_1..P_j) where the ears are paths with interiors
     DISJOINT from C, pairwise vertex-disjoint, at most one ear-foot
     per vertex, and
        |C'| = |C| + sum d_i - (removed C-edges),
     with #ears = #kept segments = #removed segments.  Every circuit
     meeting C is ONE generalized composite step from C.  Machine:
     every meeting circuit pair of K4, Petersen, Mobius-Kantor,
     Heawood - full structural verification and length formula.
 (C) ISLAND CENSUS.  Per graph and level t: the v_2 = t circuits
     ("the island"); how many are WEAK-shallow (no ear shift of
     valuation exactly t - the necessary condition 432.1 imposes on
     every island member of a residual world) and how many are
     STRONG-shallow (all shifts of valuation < t, the 431 census
     notion).
 (D) THE PAIR BATTERY.  For every strong-shallow cycle found: sweep
     all vertex-disjoint ear pairs and all four arc choices; a
     composite that is a single cycle with v_2(length) > t is a
     PROGRESS KILL (432.1, verbatim).  If pairs fail, sweep odd
     triples.  Question: does any censused shallow cycle survive
     even one generation of the composite walk?
 (E) MIXING DATA.  For a shallow base cycle of Mobius-Kantor: the
     distribution of |C'| mod 16 over all circuits meeting C - the
     first measurement of the exchange walk's character mixing.
"""
import itertools, sys, time
sys.setrecursionlimit(100000)

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

def edge_set(cyc):
    return {frozenset((a, b)) for a, b in zip(cyc, cyc[1:] + cyc[:1])}

def v2(d):
    if d == 0: return 99
    d = abs(d); r = 0
    while d % 2 == 0: d //= 2; r += 1
    return r

def ears_of(adj, C):
    Cset = set(C); L = len(C); res = []
    for i in range(L):
        for j in range(i + 1, L):
            u, v = C[i], C[j]
            stack = [(u, [u])]
            while stack:
                w, path = stack.pop()
                for nb in adj[w]:
                    if nb == v:
                        if len(path) == 1 and (j - i == 1 or
                                               (i == 0 and j == L - 1)):
                            continue
                        res.append((i, j, path + [v]))
                    elif nb not in Cset and nb not in path:
                        stack.append((nb, path + [nb]))
    return res

def decompose(C, Cp):
    """One-step decomposition of Cp relative to C.  Returns dict of
    checks + ear data, or 'disjoint'."""
    Vc = set(C)
    if not (set(Cp) & Vc): return 'disjoint'
    Ec = edge_set(C)
    m = len(Cp)
    ed = [frozenset((Cp[i], Cp[(i + 1) % m])) for i in range(m)]
    inC = [e in Ec for e in ed]
    # (A)(ii): every shared vertex uses a C-edge on Cp
    okA = True
    for i, v in enumerate(Cp):
        if v in Vc:
            e_prev = ed[(i - 1) % m]; e_next = ed[i]
            if not (e_prev in Ec or e_next in Ec): okA = False
    if all(inC):
        return {'okA': okA, 'ears': [], 'kept': len(ed), 'trivial': True,
                'formula': len(Cp) == len(C), 'structure': True}
    # rotate so position 0 starts a non-C run boundary
    s = 0
    while s < m and not (inC[(s - 1) % m] and not inC[s]): s += 1
    if s == m:
        # Cp has no C-edges at all - impossible when vertices meet
        return {'okA': okA, 'ears': None, 'structure': False,
                'formula': False, 'trivial': False}
    idx = [(s + i) % m for i in range(m)]
    runs = []          # (is_ear, [vertex list])
    cur = [Cp[idx[0]]]; cur_ear = True
    for i in range(m):
        e = ed[idx[i]]
        nxt = Cp[idx[(i + 1) % m]]
        if (e in Ec) == (not cur_ear):
            cur.append(nxt)
        else:
            runs.append((cur_ear, cur))
            cur_ear = not cur_ear
            cur = [cur[-1], nxt] if False else [Cp[idx[i]], nxt]
    runs.append((cur_ear, cur))
    # merge: runs alternate ear / kept
    ears = [r[1] for r in runs if r[0]]
    kepts = [r[1] for r in runs if not r[0]]
    # structural checks
    ok = okA
    feet = []
    for P in ears:
        if len(P) < 2: ok = False; continue
        if not (P[0] in Vc and P[-1] in Vc): ok = False
        for w in P[1:-1]:
            if w in Vc: ok = False
        feet.append(P[0]); feet.append(P[-1])
    if len(set(feet)) != len(feet): ok = False
    # pairwise vertex-disjoint ears
    seen = set()
    for P in ears:
        for w in P:
            if w in seen: ok = False
            seen.add(w)
    kept_edges = sum(inC)
    removed = len(C) - kept_edges
    dsum = sum(len(P) - 1 for P in ears)
    formula = (len(Cp) == len(C) + dsum - removed)
    counts = (len(ears) == len(kepts))
    return {'okA': okA, 'ears': ears, 'structure': ok and counts,
            'formula': formula, 'trivial': False}

GRAPHS = [('K4', k4(), 2), ('Petersen', petersen(), 3),
          ('Mobius-Kantor', mobius_kantor(), 2),
          ('Heawood', heawood(), 2)]

print('=== (A)+(B) ONE-STEP THEOREM, full verification ===', flush=True)
CENSUS = {}
for name, nE, tlv in GRAPHS:
    n, E, adj = build(nE)
    cycles = all_cycles(adj, n)
    CENSUS[name] = (n, E, adj, cycles)
    meet = disj = 0
    okA = okS = okF = 0
    t0 = time.time()
    for C in cycles:
        for Cp in cycles:
            if C == Cp: continue
            r = decompose(list(C), list(Cp))
            if r == 'disjoint':
                disj += 1; continue
            meet += 1
            okA += r['okA']; okS += r['structure']; okF += r['formula']
    print('%-14s %4d circuits: meeting pairs %6d (all: shared-edge %d,'
          ' structure %d, length formula %d) | vertex-disjoint pairs %d'
          '  [%.0fs]'
          % (name, len(cycles), meet, okA, okS, okF, disj,
             time.time() - t0), flush=True)

print(); print('=== (C) ISLAND CENSUS ===', flush=True)
SHALLOW = {}
for name, nE, tlv in GRAPHS:
    n, E, adj, cycles = CENSUS[name][0], CENSUS[name][1], \
        CENSUS[name][2], CENSUS[name][3]
    for t in (2, 3):
        isl = [C for C in cycles if v2(len(C)) == t]
        if not isl: continue
        weak = strong = 0
        stronglist = []
        for C in isl:
            shifts = []
            L = len(C)
            for (i, j, P) in ears_of(adj, C):
                d = len(P) - 1; a1 = j - i; a2 = L - a1
                shifts.append(d - a1); shifts.append(d - a2)
            if all(v2(s) != t for s in shifts): weak += 1
            if all(v2(s) < t for s in shifts):
                strong += 1; stronglist.append(C)
        print('%-14s t=%d: island size %3d | weak-shallow %d |'
              ' strong-shallow %d'
              % (name, t, len(isl), weak, strong), flush=True)
        if stronglist:
            SHALLOW[(name, t)] = stronglist

print(); print('=== (D) THE PAIR BATTERY on every strong-shallow cycle ===',
               flush=True)
def surgery(C, earA, earB, adj):
    """Composite C (+) earA-cycle (+) earB-cycle for all arc choices;
    yield (single?, length, v2)."""
    L = len(C)
    Ec = edge_set(C)
    out = []
    def arc_edges(i, j, which):
        # arc from position i to j forward (which=0) or complement
        eds = set()
        if which == 0:
            p = i
            while p != j:
                eds.add(frozenset((C[p], C[(p + 1) % L])))
                p = (p + 1) % L
        else:
            p = j
            while p != i:
                eds.add(frozenset((C[p], C[(p + 1) % L])))
                p = (p + 1) % L
        return eds
    def path_edges(P):
        return {frozenset((a, b)) for a, b in zip(P, P[1:])}
    (iA, jA, PA), (iB, jB, PB) = earA, earB
    for wA in (0, 1):
        for wB in (0, 1):
            zA = arc_edges(iA, jA, wA) | path_edges(PA)
            zB = arc_edges(iB, jB, wB) | path_edges(PB)
            comp = (Ec ^ zA) ^ zB
            # single-cycle test
            deg = {}
            for e in comp:
                a, b = tuple(e)
                deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
            if any(d != 2 for d in deg.values()): continue
            verts = list(deg)
            seen = {verts[0]}; st = [verts[0]]
            adj2 = {}
            for e in comp:
                a, b = tuple(e)
                adj2.setdefault(a, []).append(b)
                adj2.setdefault(b, []).append(a)
            while st:
                x = st.pop()
                for y in adj2[x]:
                    if y not in seen: seen.add(y); st.append(y)
            if len(seen) != len(verts): continue
            out.append(len(comp))
    return out

for (name, t), lst in sorted(SHALLOW.items()):
    n, E, adj, cycles = CENSUS[name]
    for C in lst:
        L = len(C)
        ears = ears_of(adj, list(C))
        # vertex-disjoint ear pairs
        killed = None
        npairs = 0
        for a in range(len(ears)):
            for b in range(a + 1, len(ears)):
                (i1, j1, P1), (i2, j2, P2) = ears[a], ears[b]
                if set(P1) & set(P2): continue
                npairs += 1
                for ln in surgery(list(C), ears[a], ears[b], adj):
                    if v2(ln) > t:
                        killed = (ears[a], ears[b], ln); break
                if killed: break
            if killed: break
        status = ('KILLED by pair: composite length %d (v2=%d)'
                  % (killed[2], v2(killed[2]))) if killed else \
                 'SURVIVES all %d pair composites' % npairs
        print('%-14s t=%d cycle %s: %s'
              % (name, t, ''.join(str(x) for x in C) if len(C) < 13
                 else ('len %d' % L), status), flush=True)

print(); print('=== (E) MIXING DATA (Mobius-Kantor, shallow base) ===',
               flush=True)
n, E, adj, cycles = CENSUS['Mobius-Kantor']
base = SHALLOW.get(('Mobius-Kantor', 2), [None])[0]
if base:
    hist = {}
    meet = 0
    for Cp in cycles:
        if Cp == base: continue
        r = decompose(list(base), list(Cp))
        if r == 'disjoint': continue
        meet += 1
        hist[len(Cp) % 16] = hist.get(len(Cp) % 16, 0) + 1
    print('base = shallow 12-cycle; circuits meeting base: %d of %d'
          % (meet, len(cycles) - 1))
    print('|C\'| mod 16 histogram: %s'
          % {k: hist[k] for k in sorted(hist)})
    print('(0 mod 16 would be v2 >= 4; 8 mod 16 is v2 = 3: both are'
          ' PROGRESS classes above the base t = 2)')
