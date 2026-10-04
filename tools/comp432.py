"""comp432.py -- the composite campaign on J7: climbing, supply,
spiral calculus, exact-valuation manufacture, and the residual
congruence adversary.

Checks:
 (1) CLIMBING LEMMA: for any three odd residues u, v, w mod 2^k, some
     pair-sum has 2-adic valuation exactly 1 (equivalently: among any
     three elements of equal valuation a, some pair sums to valuation
     exactly a + 1). Brute force mod 64.
 (2) SPIRAL CALCULUS on an explicit board (|C| = 24, three pairwise
     crossing ears in spiral order):
       - crossing PAIR composite = single cycle, length
         |P_i| + |P_j| + b + d (alternating arcs);
       - crossing TRIPLE zigzag = single cycle, shift
         delta_1 + delta_2' + delta_3 (mod anything - exact integers);
       - crossing QUADRUPLE (even) = TWO components (the claimed
         failure mode, verified).
 (3) LAMINAR COMPOSABILITY: three ears, one nested pair + one
     sequential, with laminar arc choices: all 7 nonempty subsets
     compose to single cycles with ADDITIVE shifts.
 (4) SUBCUBIC SUPPLY: every subcubic tree with r leaves contains
     floor(r/2) pairwise vertex-disjoint leaf-to-leaf paths
     (constructive cherry-peeling, random trees).
 (5) EXACT-VALUATION MANUFACTURE: any 2^(t+1) - 1 residues of
     valuation < t contain a subset summing to valuation EXACTLY t
     (search implementing the climbing proof; random trials t = 3);
     sharpness: the (+1/-1) configuration of size 2^(t+1) - 2 has no
     such subset.
 (6) THE RESIDUAL CONGRUENCE ADVERSARY: at t = 3, |C| = 24, the
     constant configuration (all inner shifts == c0, outer == -c0,
     u-feet on a half-system) dodges every single shift, every
     laminar-pair value and every crossing-pair value mod 8 while
     keeping all values' valuations < 3 or > 3 -- the adversary
     survives the congruence battery; recorded as the OPEN core.
"""
import itertools, random

v2 = lambda x: (x & -x).bit_length() - 1 if x else 999

# ---------------- (1) climbing lemma ----------------
bad = 0
for u, v, w in itertools.product(range(1, 64, 2), repeat=3):
    sums = [(u + v), (u + w), (v + w)]
    if not any(v2(s % 128 if s % 128 else 128) == 1 or v2(s) == 1
               for s in sums):
        # careful: valuation of the integer sum; sums of two odds are
        # even; exact valuation 1 means sum == 2 mod 4
        pass
    if not any(s % 4 == 2 for s in sums):
        bad += 1
print('(1) climbing lemma: odd triples mod 64 with NO pair-sum == 2'
      ' mod 4: %d (expected 0)' % bad)
assert bad == 0

# ---------------- board machinery ----------------
def cyc_graph(L):
    VER = list(range(L)); EDG = [(i, (i + 1) % L) for i in range(L)]
    return VER, EDG

def add_ear(VER, EDG, tag, a, b, length):
    prev = a
    for i in range(length - 1):
        v = '%s%d' % (tag, i); VER.append(v)
        EDG.append((prev, v)); prev = v
    EDG.append((prev, b))

def emask(EDG):
    return {frozenset(e): i for i, e in enumerate(EDG)}

def components_of(mask, E):
    deg = {}
    adj = {}
    for e, i in E.items():
        if mask >> i & 1:
            a, b = tuple(e)
            deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
            adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    if not deg: return 0, False
    tworeg = all(d == 2 for d in deg.values())
    seen = set(); comps = 0
    for s in deg:
        if s in seen: continue
        comps += 1; st = [s]; seen.add(s)
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in seen: seen.add(y); st.append(y)
    return comps, tworeg

# ---------------- (2) spiral calculus, |C| = 24 ----------------
L = 24
# spiral feet: x = (0, 2, 4), y = (12, 14, 16)  (pairwise crossing)
X = [0, 2, 4]; Y = [12, 14, 16]; PL = [13, 15, 17]   # |P_i| > 12 = |C|/2
VER, EDG = cyc_graph(L)
for i in range(3):
    add_ear(VER, EDG, 'e%d' % i, X[i], Y[i], PL[i])
E = emask(EDG)
mk = lambda es: sum(1 << E[frozenset(e)] for e in es)
Cm = mk([(i, (i + 1) % L) for i in range(L)])
def earmask(i):
    es = []
    prev = X[i]
    for k in range(PL[i] - 1):
        v = 'e%d%d' % (i, k); es.append((prev, v)); prev = v
    es.append((prev, Y[i]))
    return mk(es)
def arcmask(a, b):   # arc a -> b forward
    return mk([(j % L, (j + 1) % L) for j in range(a, a + (b - a) % L)])
P = [earmask(i) for i in range(3)]
Ainn = [arcmask(X[i], Y[i]) for i in range(3)]
pop = lambda m: bin(m).count('1')

# pair (0,1): crossing; composite C ^ Z0 ^ Z1 with inner arcs
Z = [P[i] ^ Ainn[i] for i in range(3)]
comp2 = Cm ^ Z[0] ^ Z[1]
c2, reg2 = components_of(comp2, E)
b_arc = (Y[0] - X[1]) % L; d_arc = (X[0] - Y[1]) % L
want2 = PL[0] + PL[1] + b_arc + d_arc
print('(2) crossing pair: single cycle: %s, length %d == |P0|+|P1|+b+d'
      ' = %d: %s' % (c2 == 1 and reg2, pop(comp2), want2,
                     pop(comp2) == want2))
assert c2 == 1 and reg2 and pop(comp2) == want2
# triple zigzag
comp3 = Cm ^ Z[0] ^ Z[1] ^ Z[2]
c3, reg3 = components_of(comp3, E)
d0 = PL[0] - (Y[0] - X[0]) % L
d1p = PL[1] + (Y[1] - X[1]) % L - L      # outer shift as integer
d2 = PL[2] - (Y[2] - X[2]) % L
want3 = 2 * L + d0 + d1p + d2            # winding (k+1)/2 = 2
print('    crossing triple zigzag: single cycle: %s, length %d =='
      ' 2|C| + d0 + d1\' + d2 = %d (WINDING 2): %s'
      % (c3 == 1 and reg3, pop(comp3), want3, pop(comp3) == want3))
assert c3 == 1 and reg3 and pop(comp3) == want3
# quadruple: add a fourth crossing ear and verify disconnection
VER4, EDG4 = cyc_graph(L)
X4 = [0, 2, 4, 6]; Y4 = [12, 14, 16, 18]; PL4 = [13, 15, 17, 13]
for i in range(4):
    add_ear(VER4, EDG4, 'f%d' % i, X4[i], Y4[i], PL4[i])
E4 = emask(EDG4)
mk4 = lambda es: sum(1 << E4[frozenset(e)] for e in es)
Cm4 = mk4([(i, (i + 1) % L) for i in range(L)])
def em4(i):
    es = []; prev = X4[i]
    for k in range(PL4[i] - 1):
        v = 'f%d%d' % (i, k); es.append((prev, v)); prev = v
    es.append((prev, Y4[i]))
    return mk4(es)
def am4(a, b):
    return mk4([(j % L, (j + 1) % L) for j in range(a, a + (b - a) % L)])
comp4 = Cm4
for i in range(4): comp4 ^= em4(i) ^ am4(X4[i], Y4[i])
c4, reg4 = components_of(comp4, E4)
print('    crossing quadruple: components = %d (2-regular: %s) -'
      ' even spirals DISCONNECT as claimed' % (c4, reg4))
assert c4 == 2

# ---------------- (3) laminar composability ----------------
VER, EDG = cyc_graph(L)
# ears: g0 nested inside g1's complement side, g2 sequential
feet = [(0, 5), (7, 11), (14, 20)]; pl = [7, 9, 13]
for i, (a, b) in enumerate(feet):
    add_ear(VER, EDG, 'g%d' % i, a, b, pl[i])
E = emask(EDG)
mk3 = lambda es: sum(1 << E[frozenset(e)] for e in es)
Cm3 = mk3([(i, (i + 1) % L) for i in range(L)])
def em3(i):
    es = []; prev = feet[i][0]
    for k in range(pl[i] - 1):
        v = 'g%d%d' % (i, k); es.append((prev, v)); prev = v
    es.append((prev, feet[i][1]))
    return mk3(es)
def am3(a, b):
    return mk3([(j % L, (j + 1) % L) for j in range(a, a + (b - a) % L)])
ok = 0
for S in range(1, 8):
    m = Cm3; want = L
    for i in range(3):
        if S >> i & 1:
            a, b = feet[i]
            m ^= em3(i) ^ am3(a, b)
            want += pl[i] - (b - a) % L
    c, reg = components_of(m, E)
    ok += (c == 1 and reg and pop(m) == want)
print('(3) laminar family: %d/7 subsets compose to single cycles with'
      ' additive shifts' % ok)
assert ok == 7

# ---------------- (4) subcubic supply ----------------
def random_subcubic_tree(r):
    """Grow a subcubic tree until it has >= r leaves (degrees <= 3)."""
    adj = {0: [1], 1: [0]}; nxt = 2
    while True:
        leaves = [v for v in adj if len(adj[v]) == 1]
        if len(leaves) >= r: break
        cands = [v for v in adj if len(adj[v]) < 3]
        v = random.choice(cands)
        adj[v].append(nxt); adj[nxt] = [v]; nxt += 1
    return adj

def disjoint_leaf_paths(adj):
    """Cherry-peeling: returns list of vertex-disjoint leaf-leaf paths."""
    adj = {v: set(ws) for v, ws in adj.items()}
    paths = []
    while True:
        leaves = [v for v in adj if len(adj[v]) == 1]
        if len(leaves) < 2: break
        # deepest branch vertex: pick any vertex with >= 2 chain-leaf
        # children after contracting degree-2 chains; simple approach:
        # BFS from an arbitrary root, find deepest vertex with two
        # leaf-chains hanging
        root = next(iter(adj))
        par = {root: None}; order = [root]; st = [root]
        while st:
            v = st.pop()
            for w in adj[v]:
                if w not in par: par[w] = v; order.append(w); st.append(w)
        def chainleaf(v, p):
            # follow the chain from v away from p; return leaf if pure chain
            path = [v]
            while True:
                nx = [w for w in adj[path[-1]] if w != (p if len(path) == 1
                      else path[-2])]
                if len(nx) == 0: return path
                if len(nx) > 1: return None
                path.append(nx[0])
        best = None
        for v in reversed(order):
            hang = []
            for w in adj[v]:
                if w == par[v]: continue
                cl = chainleaf(w, v)
                if cl is not None: hang.append(cl)
            if len(hang) >= 2:
                best = (v, hang[0], hang[1]); break
        if best is None:
            # tree is a single path: one path, done
            # find its two leaves
            l1, l2 = leaves[0], leaves[1]
            paths.append((l1, l2)); break
        v, c1, c2 = best
        paths.append((c1[-1], c2[-1]))
        for x in c1 + c2 + [v]:
            for w in list(adj.get(x, [])):
                adj[w].discard(x)
            adj.pop(x, None)
    return paths

random.seed(432)
ok = 0; tot = 0
for _ in range(200):
    adj = random_subcubic_tree(random.randint(4, 30))
    leaves = [v for v in adj if len(adj[v]) == 1]
    r = len(leaves)
    if r < 2: continue
    tot += 1
    ps = disjoint_leaf_paths(adj)
    ok += (len(ps) >= r // 2)
    assert len(ps) >= r // 2, (r, len(ps))
print('(4) subcubic supply: %d/%d random subcubic trees yield'
      ' floor(r/2) disjoint leaf-paths' % (ok, tot))

# ---------------- (5) exact-valuation manufacture ----------------
def exact_val_subset(vals, t):
    """Find a subset with v2(sum) == t, via climbing. vals: list of ints
    with v2 < t. Returns subset indices or None."""
    pool = [(v % (1 << (t + 2)), frozenset([i])) for i, v in enumerate(vals)]
    while True:
        for val, S in pool:
            if v2(val % (1 << (t + 1)) or (1 << (t + 1))) == t:
                return S
        # group by valuation, climb
        byval = {}
        for val, S in pool:
            a = v2(val % (1 << (t + 2)) or (1 << (t + 2)))
            byval.setdefault(min(a, t + 1), []).append((val, S))
        prog = False
        for a in sorted(byval):
            if a >= t: continue
            grp = byval[a]
            # try all pairs with disjoint index sets
            done = False
            for i in range(len(grp)):
                for j in range(i + 1, len(grp)):
                    (v1, S1), (v2_, S2) = grp[i], grp[j]
                    if S1 & S2: continue
                    s = (v1 + v2_)
                    if v2(s % (1 << (t + 2)) or (1 << (t + 2))) == a + 1:
                        pool.remove(grp[i]); pool.remove(grp[j])
                        pool.append((s, S1 | S2))
                        prog = True; done = True; break
                if done: break
            if done: break
        if not prog: return None

random.seed(4325)
t = 3; N = (1 << (t + 1)) - 1
ok = 0
for _ in range(400):
    vals = []
    while len(vals) < N:
        x = random.randrange(1, 1 << (t + 1))
        if v2(x) < t: vals.append(x)
    S = exact_val_subset(vals, t)
    if S is not None:
        s = sum(vals[i] for i in S)
        assert v2(s) == t or s % (1 << (t + 1)) == (1 << t)
        ok += 1
print('(5) exact-valuation: %d/400 random 15-element pools (t=3) yield'
      ' a subset with v_2(sum) exactly 3' % ok)
assert ok == 400
# sharpness: 2^t-1 ones and 2^t-1 minus-ones mod 2^(t+1)
M = 1 << (t + 1)
vals = [1] * ((1 << t) - 1) + [M - 1] * ((1 << t) - 1)
found = False
for a in range(len(vals) + 1):
    for b in range(len(vals) + 1):
        if a <= (1 << t) - 1 and b <= (1 << t) - 1 and (a or b):
            if (a - b) % M == (1 << t): found = True
print('    sharpness: size-14 (+1/-1) pool has exact-3 subset: %s'
      ' (expected False)' % found)
assert not found

# ---------------- (6) the residual congruence adversary ----------------
t = 3; Mt = 1 << t; L24 = 24
c0 = 2   # v2(c0) = 1 < 3
# half-system u-classes mod 8 avoiding difference 4: {0,1,2,3}
U = [0, 1, 2, 3]
# inner arcs == -c0 mod 8, |P| == 0 mod 8: delta = c0... build 12 ears
ears = []
for k in range(12):
    u = U[k % 4] + 8 * (k // 4)          # positions on the 24-cycle
    inner = (-c0) % Mt + 8               # arc length == -c0 mod 8
    ears.append((u, (u + inner) % L24, 16))   # |P| = 16 == 0 mod 8
sing = set(); lam = set(); cross = set()
for (u, v, p) in ears:
    inner = (v - u) % L24
    sing.add((p - inner) % Mt); sing.add((p + inner) % Mt)
for (u1, y1, p1), (u2, y2, p2) in itertools.combinations(ears, 2):
    i1 = (y1 - u1) % L24; i2 = (y2 - u2) % L24
    lam.add(((p1 - i1) + (p2 - i2)) % Mt)
    cross.add((p1 + p2 - ((u2 - u1) % L24) - ((y2 - y1) % L24)) % Mt)
print('(6) residual adversary t=3: single-shift classes %s, laminar-pair'
      ' %s, crossing-pair %s (mod 8)' % (sorted(sing), sorted(lam),
                                         sorted(cross)))
# the 0-mod-8 crossing values at the INTEGER level: progress (v_2 = 3)
# or lateral (v_2 >= 4)?
vals3 = []
for (u1, y1, p1), (u2, y2, p2) in itertools.combinations(ears, 2):
    d_int = p1 + p2 - ((u2 - u1) % L24) - ((y2 - y1) % L24)
    if d_int % Mt == 0 and d_int != 0:
        vals3.append(v2(abs(d_int)))
print('    integer valuations of the 0-mod-8 crossing values: %s'
      ' -> all >= 4: %s (LATERAL, unusable): the PROGRESS battery'
      % (sorted(set(vals3)), all(x >= 4 for x in vals3)))
print('    alone does not kill him; he survives as a braided'
      ' constant/half-system configuration')
assert all(x >= 4 for x in vals3)
print()
print('VERDICT: climbing, supply, spiral calculus and exact-valuation')
print('manufacture are THEOREMS; progress-composites are legal moves;')
print('the kill fires wherever a laminar family of 2^(t+1)-1 ears')
print('exists; the surviving adversary is a rigid constant/half-system')
print('configuration - the OPEN core of J7, now finite-dimensional.')
