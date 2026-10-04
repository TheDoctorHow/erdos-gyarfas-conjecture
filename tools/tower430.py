"""tower430.py -- naming the machine: the Brown form, divisible codes,
and the polarization tower.

Machine checks for section 430:
 (1) INCLUSION-EXCLUSION: |x1 ^ ... ^ xk| = sum_{S nonempty} (-2)^(|S|-1)
     |intersection over S|  (exact, all k) -- the master identity behind
     every level of the tower.
 (2) DEGREE-t STRUCTURE: the k-th finite difference of q_t(x) = |x| mod 2^t
     equals (-2)^(k-1) |x1 n ... n xk| mod 2^t; hence differences of order
     > t vanish: q_t is a DEGREE-t form.  t=1 linear (bipartite parity),
     t=2 quadratic (the Brown form), t=3 cubic (the phase boundary).
 (3) DIVISIBILITY CASCADE on codes: a triply even code (weights == 0 mod 8)
     has pairwise intersections == 0 mod 4 and triple intersections
     == 0 mod 2; Ward's crude binary corollary k <= (e+1) n / 2^e checked.
 (4) WEB POLARIZATION DESCENT: on random bridges, even connecting paths
     aligned mod 8 have |P ^ P'| == 2(l0 - |P n P'|) mod 8: level-t length
     rigidity IS level-(t-1) rigidity of the intersection function.
 (5) BROWN RIGIDITY vs CUBIC WILDNESS: Gauss sums over random subcodes of
     cycle spaces: at t=2, |G|^2 in {0} u {4^d 2^rho} (Brown's theorem,
     with d = dim radical, rho = rank); at t=3 (zeta_8) the magnitude
     census is DIVERSE -- no Brown law survives.
 (6) RADICAL = BICYCLE SPACE: on graph cycle spaces the radical of the
     intersection pairing equals C n C-perp (every radical element lies in
     the cut space), and the Gauss law holds with d = bicycle dimension.
"""
import random, itertools, cmath

pop = lambda m: bin(m).count('1')

# ---------------- (1) inclusion-exclusion ----------------
random.seed(430)
UNIV = 40
def randset():
    return sum(1 << i for i in range(UNIV) if random.random() < 0.4)

ok = 0; tot = 0
for k in (2, 3, 4, 5):
    for _ in range(60):
        xs = [randset() for _ in range(k)]
        xor = 0
        for x in xs: xor ^= x
        s = 0
        for r in range(1, k + 1):
            for S in itertools.combinations(range(k), r):
                cap = (1 << UNIV) - 1
                for i in S: cap &= xs[i]
                s += (-2) ** (r - 1) * pop(cap)
        tot += 1; ok += (pop(xor) == s)
        assert pop(xor) == s
print('(1) inclusion-exclusion |x1^...^xk| = sum (-2)^(|S|-1)|cap S|:'
      ' %d/%d exact' % (ok, tot))

# ---------------- (2) degree-t structure of q_t ----------------
def diff(f, xs):
    # k-th finite difference of f at 0 along x1..xk (integer values)
    k = len(xs); s = 0
    for r in range(k + 1):
        for S in itertools.combinations(range(k), r):
            v = 0
            for i in S: v ^= xs[i]
            s += (-1) ** (k - r) * f(v)
    return s

t = 3; MOD = 1 << t
ok = 0; tot = 0
for k in (1, 2, 3, 4, 5):
    for _ in range(40):
        xs = [randset() for _ in range(k)]
        cap = (1 << UNIV) - 1
        for x in xs: cap &= x
        d = diff(lambda v: pop(v), xs) % MOD
        want = ((-2) ** (k - 1) * pop(cap)) % MOD if k >= 1 else 0
        tot += 1; ok += (d == want)
        assert d == want
        if k > t:
            assert d % MOD == ((-2) ** (k - 1) * pop(cap)) % MOD == 0 \
                if (1 << (k - 1)) % MOD == 0 else True
print('(2) k-th difference of q_%d == (-2)^(k-1)|cap| mod %d: %d/%d;'
      ' order > %d differences vanish => q_%d is a DEGREE-%d form'
      % (t, MOD, ok, tot, t, t, t))

# ---------------- (3) divisibility cascade on codes ----------------
def greedy_doubly_even(n, tries=300):
    gens = []
    for _ in range(tries):
        w = 0
        idx = random.sample(range(n), random.choice([4, 8]))
        for i in idx: w |= 1 << i
        if pop(w) % 4: continue
        if any(pop(w & g) % 2 for g in gens): continue
        gens.append(w)
    return gens

def span(gens):
    out = {0}
    for g in gens:
        out |= {x ^ g for x in out}
    return out

casc_ok = True
for trial in range(6):
    n0 = 24
    gens = greedy_doubly_even(n0)
    C0 = span(gens)
    assert all(pop(x) % 4 == 0 for x in C0)
    # replicate: (c, c) -- triply even of length 2 n0
    C1 = {x | (x << n0) for x in C0}
    n1, k1 = 2 * n0, len(gens)
    assert all(pop(x) % 8 == 0 for x in C1)
    words = [x for x in C1 if x][:14]
    for a, b in itertools.combinations(words, 2):
        casc_ok &= (pop(a & b) % 4 == 0)
    for a, b, c in itertools.combinations(words[:9], 3):
        casc_ok &= (pop(a & b & c) % 2 == 0)
    # Ward crude corollary, e = 3 (weights all == 0 mod 8)
    assert k1 <= 4 * n1 // 8, 'Ward corollary violated'
assert casc_ok
print('(3) triply even codes (replicated doubly even, n=48): weights == 0'
      ' mod 8, pairwise n == 0 mod 4, triple n == 0 mod 2: CASCADE holds;'
      ' k <= (e+1)n/2^e checked')

# ---------------- (4) web polarization descent on bridges ----------------
def random_bridge(n):
    for _ in range(4000):
        stubs = [v for v in range(n) for _ in range(3)]
        random.shuffle(stubs)
        E = set(); good = True
        for i in range(0, len(stubs), 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b or frozenset((a, b)) in E: good = False; break
            E.add(frozenset((a, b)))
        if not good: continue
        adj = {v: [] for v in range(n)}
        for e in E:
            a, b = tuple(e); adj[a].append(b); adj[b].append(a)
        g = 99
        for s in range(n):
            dist = {s: 0}; par = {s: None}; q = [s]
            while q:
                v = q.pop(0)
                for w in adj[v]:
                    if w not in dist:
                        dist[w] = dist[v] + 1; par[w] = v; q.append(w)
                    elif par[v] != w: g = min(g, dist[v] + dist[w] + 1)
        if g < 5: continue
        El = [tuple(e) for e in E]; random.shuffle(El)
        picked = []
        for e in El:
            if all(v not in f for f in picked for v in e): picked.append(e)
            if len(picked) == 3: break
        if len(picked) < 3: continue
        VER = list(range(n)) + ['t0', 't1', 't2']
        EDG = [tuple(e) for e in E
               if tuple(sorted(e)) not in [tuple(sorted(x)) for x in picked]]
        for kk, (a, b) in enumerate(picked):
            EDG += [(a, 't%d' % kk), ('t%d' % kk, b)]
        return VER, EDG
    return None

def machinery(VER, EDG):
    E = {frozenset(e): i for i, e in enumerate(EDG)}
    adj = {v: [] for v in VER}
    for u, w in EDG: adj[u].append(w); adj[w].append(u)
    bit = lambda u, w: 1 << E[frozenset((u, w))]
    return E, adj, bit

def all_paths(adj, bit, src, dst, cap=4000):
    out = []
    def dfs(v, mask, seen):
        if len(out) >= cap: return
        if v == dst: out.append(mask); return
        for w in adj[v]:
            if w not in seen: dfs(w, mask | bit(v, w), seen | {w})
    dfs(src, 0, {src})
    return out

pairs8 = 0; desc_ok = 0; iota_classes = {}
for n in (12, 14, 16):
    res = random_bridge(n)
    if not res: continue
    VER, EDG = res
    E, adj, bit = machinery(VER, EDG)
    for (s, tt) in [('t0', 't1'), ('t0', 't2'), ('t1', 't2')]:
        ev = [m for m in all_paths(adj, bit, s, tt) if pop(m) % 2 == 0]
        random.shuffle(ev); ev = ev[:60]
        for i in range(len(ev)):
            for j in range(i + 1, min(i + 12, len(ev))):
                P, Q = ev[i], ev[j]
                if pop(P) % 8 != pop(Q) % 8: continue
                pairs8 += 1
                l0 = pop(P); iota = pop(P & Q)
                lhs = pop(P ^ Q) % 8
                rhs = (2 * (l0 - iota)) % 8
                desc_ok += (lhs == rhs)
                assert lhs == rhs
                iota_classes[iota % 4] = iota_classes.get(iota % 4, 0) + 1
print('(4) web descent |P^Q| == 2(l0 - |PnQ|) mod 8 on mod-8-aligned pairs:'
      ' %d/%d;  iota mod 4 census: %s' % (desc_ok, pairs8,
      dict(sorted(iota_classes.items()))))
print('    => pinning |P^Q| mod 8 is EXACTLY pinning |PnQ| mod 4:'
      ' level-3 length rigidity = level-2 intersection rigidity')

# ---------------- (5) Brown rigidity vs cubic wildness ----------------
def rank_radical(basis):
    k = len(basis)
    gram = [[pop(basis[i] & basis[j]) % 2 for j in range(k)]
            for i in range(k)]
    # F2 Gaussian elimination for rank
    M = [row[:] for row in gram]; r = 0
    for c in range(k):
        piv = next((i for i in range(r, k) if M[i][c]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(k):
            if i != r and M[i][c]:
                M[i] = [(a + b) % 2 for a, b in zip(M[i], M[r])]
        r += 1
    return r, k - r

def basis_of(words):
    basis = []
    for w in words:
        x = w
        for b in basis: x = min(x, x ^ b)
        red = w
        for b in sorted(basis, reverse=True):
            red = min(red, red ^ b)
        if red: basis.append(red)
    return basis

def fundcycles(VER, EDG):
    E, adj, bit = machinery(VER, EDG)
    tree = set(); seen = {VER[0]}; st = [VER[0]]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w); tree.add(frozenset((v, w))); st.append(w)
    par = {VER[0]: None}; st = [VER[0]]
    while st:
        v = st.pop()
        for w in adj[v]:
            if frozenset((v, w)) in tree and w not in par:
                par[w] = v; st.append(w)
    def tp(a, b):
        m = 0; anc = set(); x = a
        while x is not None: anc.add(x); x = par[x]
        y = b
        while y not in anc: m |= bit(y, par[y]); y = par[y]
        x = a
        while x != y: m |= bit(x, par[x]); x = par[x]
        return m
    return E, [(1 << E[e]) | tp(*tuple(e)) for e in E if e not in tree]

brown_ok = 0; brown_tot = 0; wild = set()
z8 = cmath.exp(2j * cmath.pi / 8)
for n in (10, 12, 14):
    res = random_bridge(n)
    if not res: continue
    VER, EDG = res
    E, fund = fundcycles(VER, EDG)
    for _ in range(8):
        sub = random.sample(fund, random.randint(3, min(6, len(fund))))
        basis = basis_of(sub)
        C = span(basis)
        rho, d = rank_radical(basis)
        G2 = sum(1j ** (pop(x) % 4) for x in C)
        m2 = abs(G2) ** 2
        law = {0.0} | {round((4 ** d) * (2 ** rho), 6)}
        brown_tot += 1
        okb = round(m2, 6) in law
        brown_ok += okb
        assert okb, (m2, d, rho)
        G3 = sum(z8 ** (pop(x) % 8) for x in C)
        wild.add(round(abs(G3) ** 2, 3))
print('(5) Brown law |G_2|^2 in {0, 4^d 2^rho}: %d/%d subcodes;'
      % (brown_ok, brown_tot))
print('    zeta_8 census |G_3|^2 values: %d distinct magnitudes %s'
      % (len(wild), sorted(wild)[:10]))
assert len(wild) > 4
print('    => NO Brown-type magnitude law at t = 3: the phase boundary,'
      ' now in the Weil-sum frame')

# ---------------- (6) radical = bicycle space ----------------
K4 = (list(range(4)), [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
tests = [('K4 (Tutte(-1,-1) = 4 => bicycle dim 2)', K4)]
for n in (10, 12):
    res = random_bridge(n)
    if res: tests.append(('random bridge n=%d' % n, res))
for name, (VER, EDG) in tests:
    E, fund = fundcycles(VER, EDG)
    basis = basis_of(fund)
    k = len(basis)
    # radical vectors: solve Gram x = 0 by brute force over small C
    C = span(basis)
    rad = [x for x in C if all(pop(x & b) % 2 == 0 for b in basis)]
    # cut space: generated by vertex stars
    stars = []
    for v in VER:
        m = 0
        for e, i in E.items():
            if v in e: m |= 1 << i
        stars.append(m)
    cutb = basis_of(stars)
    incut = lambda x: not [1 for b in [x] if _reduce(b, cutb)]
    def _reduce(x, bs):
        for b in sorted(bs, reverse=True):
            x = min(x, x ^ b)
        return x
    allcut = all(_reduce(x, cutb) == 0 for x in rad)
    rho, d = rank_radical(basis)
    assert len(rad) == 2 ** d
    print('(6) %s: dim C = %d, radical dim = %d, every radical elt in'
          ' cut space: %s' % (name, k, d, allcut))
    assert allcut
    if name.startswith('K4'): assert d == 2

print()
print('VERDICT: the t = 2 machine is the BROWN FORM of the cycle matroid')
print('(radical = bicycle space, Gauss law = Brown/Vertigan); the tower is')
print('the polarization tower of a DEGREE-t form; divisibility cascades one')
print('level down per polarization; Brown rigidity dies at t = 3 but the')
print('descent identity survives at every t.')
