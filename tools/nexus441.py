"""nexus441.py -- the meditation machine: every stratum t is the shadow
of ONE 2-adic quadratic form on the integral flow lattice.

 (A) WHAT EXACTLY BREAKS AT MOD 8.
     (A1) Lucas digit identity: the k-th binary digit of |Z| is
          binom(|Z|, 2^k) mod 2 = e_{2^k} of the indicator vector: the
          new form at each level is an F_2-polynomial of degree 2^k.
     (A2) Third differences: |.| mod 8 on the XOR space has nonzero
          third difference (degree >= 3); |.|^2 on the lattice has
          third difference IDENTICALLY ZERO (quadratic forever).
     (A3) Wildness table: random subspaces of F_2^12; count distinct
          Gauss magnitudes |sum zeta_{2^t}^{|x|}| at t = 2 (Brown:
          rigid) versus t = 3 (wild). The classification that closed
          t = 2 does not exist in XOR coordinates.

 (B) THE TOWER THEOREM (the lift).  Lambda = ker(oriented incidence)
     = the lattice of integral circulations.  Circuits lift to +-1
     circulations with |x|^2 = length EXACTLY (machine: full circuit
     census, four graphs).  For every s >= 1,
         q_s(x) = |x|^2 mod 2^{s+1}
     is WELL DEFINED on Lambda / 2^s Lambda and quadratic with
     polarization 2<x,y> (machine identity checks, s = 1..4).
     Archive stratum t <-> lattice level s = t - 1;  s = 1 is the
     Brown form of section 430.

 (C) INVARIANTS AT EVERY LEVEL.
     (C1) Lattice Gauss sums G_s = sum zeta_{2^{s+1}}^{q_s(x)} over
          Lambda/2^s Lambda: log2 |G_s|^2 (rigidity check).
     (C2) det Gram(fundamental cycle basis) = number of spanning
          trees (Kirchhoff: discriminant group = critical group).
     (C3) corank of Gram mod 2 = bicycle dimension (radical law,
          all four graphs) - generalizing 430's radical = bicycle.

 (D) THE REFORMULATION (Mobius-Kantor, full census).
     EAR IDENTITY:  delta = q(z) + 2<x, z>  for the ear-cycle lift z
     oriented against x on its arc (every cycle, every ear, both
     arcs, exact).
     DEEP-EAR CRITERION:  2^t | delta  <=>  q(z) even  AND
     <x,z> == -q(z)/2 (mod 2^{t-1}).
     EQUIVALENCE:  shallow_direct == shallow_lattice for EVERY cycle
     at t = 2 and t = 3; counts reported (expect 4 / 0 per 431).
"""
import cmath, itertools, math, random, sys, time
sys.setrecursionlimit(100000)
random.seed(441)

# ---------------- graphs ----------------
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
    eidx = {}
    for k, (a, b) in enumerate(E):
        eidx[(a, b)] = (k, 1); eidx[(b, a)] = (k, -1)
    return n, E, adj, eidx

def norm2(x): return sum(v * v for v in x)
def dot(x, y): return sum(a * b for a, b in zip(x, y))
def v2(d):
    if d == 0: return 99
    d = abs(d); r = 0
    while d % 2 == 0: d //= 2; r += 1
    return r

# ---------------- lattice basis ----------------
def flow_basis(n, E, adj, eidx):
    """Fundamental cycle basis of the integral flow lattice."""
    par = {0: None}; depth = {0: 0}; order = [0]
    for v in order:
        for w in adj[v]:
            if w not in par:
                par[w] = v; depth[w] = depth[v] + 1; order.append(w)
    tree = set()
    for w, p in par.items():
        if p is not None:
            tree.add((min(w, p), max(w, p)))
    def tree_path(a, b):
        pa, pb = [a], [b]
        x, y = a, b
        while depth[x] > depth[y]: x = par[x]; pa.append(x)
        while depth[y] > depth[x]: y = par[y]; pb.append(y)
        while x != y:
            x = par[x]; pa.append(x)
            y = par[y]; pb.append(y)
        return pa + pb[::-1][1:][::-1] if False else pa + list(reversed(pb[:-1]))
    B = []
    m = len(E)
    for (a, b) in E:
        if (a, b) in tree: continue
        x = [0] * m
        k, s = eidx[(a, b)]; x[k] += s
        path = tree_path(b, a)          # b -> ... -> a through the tree
        for p, q in zip(path, path[1:]):
            k, s = eidx[(p, q)]; x[k] += s
        B.append(x)
    return B

def is_circulation(n, E, x):
    bal = [0] * n
    for k, (a, b) in enumerate(E):
        bal[a] -= x[k]; bal[b] += x[k]
    return all(v == 0 for v in bal)

# ---------------- integer determinant (Bareiss) ----------------
def idet(M0):
    M = [row[:] for row in M0]; n = len(M)
    if n == 0: return 1
    sign = 1; prev = 1
    for i in range(n - 1):
        if M[i][i] == 0:
            for r in range(i + 1, n):
                if M[r][i] != 0:
                    M[i], M[r] = M[r], M[i]; sign = -sign; break
            else:
                return 0
        for r in range(i + 1, n):
            for c in range(i + 1, n):
                M[r][c] = (M[r][c] * M[i][i] - M[r][i] * M[i][c]) // prev
        prev = M[i][i]
    return sign * M[n - 1][n - 1]

def spanning_trees(n, E):
    L = [[0] * n for _ in range(n)]
    for a, b in E:
        L[a][a] += 1; L[b][b] += 1
        L[a][b] -= 1; L[b][a] -= 1
    R = [row[1:] for row in L[1:]]
    return idet(R)

# ---------------- F_2 linear algebra on bitmasks ----------------
def f2rank(rows):
    piv = {}
    for v in rows:
        for p in sorted(piv, reverse=True):
            if (v >> p) & 1: v ^= piv[p]
        if v: piv[v.bit_length() - 1] = v
    return len(piv)

# ---------------- lattice Gauss sums ----------------
def gauss_sum(G, s):
    r = len(G); M = 1 << (s + 1); counts = [0] * M
    def rec(i, acc, lin):
        if i == r:
            counts[acc % M] += 1; return
        rec(i + 1, acc, lin)
        for v in range(1, 1 << s):
            acc2 = acc + v * v * G[i][i] + 2 * v * lin[i]
            lin2 = [lin[j] + v * G[i][j] for j in range(r)]
            rec(i + 1, acc2, lin2)
    rec(0, 0, [0] * r)
    z = sum(c * cmath.exp(2j * math.pi * a / M) for a, c in enumerate(counts))
    return abs(z)

# ---------------- cycle census ----------------
def all_cycles(adj, n):
    out = set()
    def dfs(start, v, visited, path):
        for w in adj[v]:
            if w == start and len(path) >= 3:
                if path[1] < path[-1]:
                    out.add(tuple(path))
            elif w > start and w not in visited:
                visited.add(w); path.append(w)
                dfs(start, w, visited, path)
                path.pop(); visited.remove(w)
    for s in range(n):
        dfs(s, s, {s}, [s])
    return sorted(out, key=lambda c: (len(c), c))

def lift(cyc, eidx, m):
    x = [0] * m
    for a, b in zip(cyc, cyc[1:] + cyc[:1]):
        k, s = eidx[(a, b)]; x[k] += s
    return x

# =====================================================================
print('=== (A) WHAT BREAKS AT MOD 8 ===', flush=True)

# (A1) Lucas digit identity
ok = 0
for _ in range(200):
    w = random.randrange(1, 5000)
    k = random.randrange(0, 7)
    if (math.comb(w, 1 << k) % 2) == ((w >> k) & 1): ok += 1
print('(A1) digit_k(|Z|) = binom(|Z|, 2^k) mod 2 (Lucas): %d/200' % ok)

# (A2) third differences
pc = lambda x: bin(x).count('1')
bad = 0
for _ in range(300):
    x, y, z = (random.getrandbits(20) for _ in range(3))
    D3 = (pc(x ^ y ^ z) - pc(x ^ y) - pc(x ^ z) - pc(y ^ z)
          + pc(x) + pc(y) + pc(z)) % 8
    if D3 != 0: bad += 1
print('(A2) XOR third difference of |.| mod 8 nonzero: %d/300 triples' % bad)
zero = 0
for _ in range(300):
    x = [random.randrange(-5, 6) for _ in range(12)]
    y = [random.randrange(-5, 6) for _ in range(12)]
    z = [random.randrange(-5, 6) for _ in range(12)]
    sxyz = [a + b + c for a, b, c in zip(x, y, z)]
    sxy = [a + b for a, b in zip(x, y)]
    sxz = [a + c for a, c in zip(x, z)]
    syz = [b + c for b, c in zip(y, z)]
    T3 = (norm2(sxyz) - norm2(sxy) - norm2(sxz) - norm2(syz)
          + norm2(x) + norm2(y) + norm2(z))
    if T3 == 0: zero += 1
print('     lattice third difference of |.|^2 == 0: %d/300 triples' % zero)

# (A3) wildness table
def rand_subspace(nbits, d):
    piv = {}
    while len(piv) < d:
        v = random.getrandbits(nbits)
        for p in sorted(piv, reverse=True):
            if (v >> p) & 1: v ^= piv[p]
        if v: piv[v.bit_length() - 1] = v
    elems = [0]
    for b in piv.values():
        elems = elems + [e ^ b for e in elems]
    return elems

mags = {2: set(), 3: set()}
for d in (4, 6, 8):
    for _ in range(40):
        elems = rand_subspace(12, d)
        ws = [pc(e) for e in elems]
        for t in (2, 3):
            M = 1 << t
            g = abs(sum(cmath.exp(2j * math.pi * (w % M) / M) for w in ws))
            mags[t].add(round(g, 4))
print('(A3) 120 random subspaces of F_2^12, dims 4/6/8:')
for t in (2, 3):
    print('     t=%d: %2d distinct Gauss magnitudes  %s' %
          (t, len(mags[t]),
           '(RIGID)' if len(mags[t]) <= 8 else '(WILD)'))

# =====================================================================
print(); print('=== (B) THE TOWER THEOREM ===', flush=True)
GRAPHS = [('K4', k4()), ('Petersen', petersen()),
          ('Mobius-Kantor', mobius_kantor()), ('Heawood', heawood())]
DATA = {}
for name, nE in GRAPHS:
    n, E, adj, eidx = build(nE)
    m = len(E)
    B = flow_basis(n, E, adj, eidx)
    r = len(B)
    assert r == m - n + 1
    assert all(is_circulation(n, E, b) for b in B)
    # well-definedness identity, s = 1..4
    idok = 0
    for _ in range(200):
        cx = [random.randrange(-3, 4) for _ in range(r)]
        cl = [random.randrange(-3, 4) for _ in range(r)]
        x = [sum(cx[i] * B[i][k] for i in range(r)) for k in range(m)]
        lam = [sum(cl[i] * B[i][k] for i in range(r)) for k in range(m)]
        good = True
        for s in range(1, 5):
            xs = [a + (1 << s) * b for a, b in zip(x, lam)]
            if (norm2(xs) - norm2(x)) % (1 << (s + 1)) != 0: good = False
        idok += good
    # circuit lifts
    cycles = all_cycles(adj, n)
    lok = 0
    for cyc in cycles:
        x = lift(cyc, eidx, m)
        if is_circulation(n, E, x) and norm2(x) == len(cyc) \
           and all(abs(v) == 1 for v in x if v != 0):
            lok += 1
    print('%-14s rank %d, %4d circuits: lift |x|^2=length %d/%d;'
          ' q_s well-defined s=1..4: %d/200'
          % (name, r, len(cycles), lok, len(cycles), idok), flush=True)
    G = [[dot(a, b) for b in B] for a in B]
    DATA[name] = (n, E, adj, eidx, m, B, r, G, cycles)

# =====================================================================
print(); print('=== (C) INVARIANTS AT EVERY LEVEL ===', flush=True)
for name, nE in GRAPHS:
    n, E, adj, eidx, m, B, r, G, cycles = DATA[name]
    dg = idet(G); st = spanning_trees(n, E)
    # bicycle dimension
    cyc_rows = []
    for b in B:
        v = 0
        for k, val in enumerate(b):
            if val % 2: v |= (1 << k)
        cyc_rows.append(v)
    cut_rows = []
    for w in range(n):
        v = 0
        for k, (a, bb) in enumerate(E):
            if a == w or bb == w: v |= (1 << k)
        cut_rows.append(v)
    rc = f2rank(cyc_rows); rq = f2rank(cut_rows)
    rall = f2rank(cyc_rows + cut_rows)
    bic = rc + rq - rall
    gram_rows = []
    for row in G:
        v = 0
        for k, val in enumerate(row):
            if val % 2: v |= (1 << k)
        gram_rows.append(v)
    corank = r - f2rank(gram_rows)
    print('%-14s det Gram = %6d  #trees = %6d  %s |'
          ' bicycle dim %d, Gram mod-2 corank %d  %s'
          % (name, dg, st, 'OK' if dg == st else 'FAIL',
             bic, corank, 'OK' if bic == corank else 'FAIL'), flush=True)

print(); print('lattice Gauss sums: log2 |G_s|^2  (levels s; stratum t = s+1)')
for name, nE in GRAPHS:
    n, E, adj, eidx, m, B, r, G, cycles = DATA[name]
    levels = (1, 2, 3) if r <= 6 else (1, 2)
    row = []
    for s in levels:
        t0 = time.time()
        g = gauss_sum(G, s)
        if g < 1e-6:
            row.append('s=%d: G=0' % s)
        else:
            l2 = 2 * math.log2(g)
            tag = 'INT' if abs(l2 - round(l2)) < 1e-6 else 'NONINT(%.3f)' % l2
            row.append('s=%d: %s %s' % (s, ('%.0f' % round(l2)), tag))
    print('%-14s %s' % (name, '   '.join(row)), flush=True)

# =====================================================================
print(); print('=== (D) THE REFORMULATION (Mobius-Kantor census) ===',
               flush=True)
n, E, adj, eidx, m, B, r, G, cycles = DATA['Mobius-Kantor']

def ears_of(C):
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
                            continue      # the C-edge itself
                        res.append((i, j, path + [v]))
                    elif nb not in Cset and nb not in path:
                        stack.append((nb, path + [nb]))
    return res

id_ok = id_all = 0
eq_ok = eq_all = 0
shallow_counts = {2: [0, 0], 3: [0, 0]}   # [direct, lattice]
mismatch = 0
for C in cycles:
    L = len(C)
    x = lift(C, eidx, m)
    ears = ears_of(C)
    pairs = []          # (delta, q, ip) for both arcs of every ear
    for (i, j, P) in ears:
        d = len(P) - 1
        zP = [0] * m
        for a, b in zip(P, P[1:]):
            k, s = eidx[(a, b)]; zP[k] += s
        arc = C[i:j + 1]
        za1 = [0] * m
        for a, b in zip(arc, arc[1:]):
            k, s = eidx[(a, b)]; za1[k] += s
        a1 = j - i; a2 = L - a1
        z1 = [p - q for p, q in zip(zP, za1)]
        za2 = [xx - q for xx, q in zip(x, za1)]
        z2 = [-(p + q) for p, q in zip(zP, za2)]
        for (z, aa) in ((z1, a1), (z2, a2)):
            delta = d - aa
            qz = norm2(z); ip = dot(x, z)
            id_all += 1
            if delta == qz + 2 * ip: id_ok += 1
            pairs.append((delta, qz, ip))
    for t in (2, 3):
        for (delta, qz, ip) in pairs:
            deep = (delta % (1 << t) == 0)
            crit = (qz % 2 == 0) and ((ip + qz // 2) % (1 << (t - 1)) == 0)
            eq_all += 1
            if deep == crit: eq_ok += 1
        sd = (v2(L) == t) and all(v2(delta) < t for (delta, _, _) in pairs)
        sl = (v2(L) == t) and not any(
            (qz % 2 == 0) and ((ip + qz // 2) % (1 << (t - 1)) == 0)
            for (_, qz, ip) in pairs)
        shallow_counts[t][0] += sd
        shallow_counts[t][1] += sl
        if sd != sl: mismatch += 1

print('EAR IDENTITY  delta = q(z) + 2<x,z>:            %d/%d exact'
      % (id_ok, id_all))
print('DEEP-EAR CRITERION (t = 2 and 3):               %d/%d match'
      % (eq_ok, eq_all))
print('shallow_direct == shallow_lattice mismatches:   %d' % mismatch)
for t in (2, 3):
    print('  t=%d shallow cycles: direct %d, lattice %d'
          % (t, shallow_counts[t][0], shallow_counts[t][1]))
print()
print('Every stratum is the shadow of one 2-adic quadratic form;')
print('shallowness = dodging one linear congruence mod 2^(t-1) per')
print('even ear-cycle - a covering problem, uniform in t.')
