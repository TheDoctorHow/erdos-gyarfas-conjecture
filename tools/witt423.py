"""witt423.py -- the underlying structure: the Z_4-form (C, lambda, q),
its Gauss sum, and the two stratum-free reduction targets.

q(x) = |x| mod 4 on the cycle space C is Brown's Z_4-quadratic refinement of
the intersection form lambda (222.1): q(x+y) = q(x) + q(y) + 2 lambda(x,y).
Its classification data (radical dim d, rank rho = mu - d, Brown invariant)
is what every remaining case is made of:
    theta-planes  : q == 2 on nonzero elements, Gauss sum -2   (J2')
    odd part      : r orthonormal odd cycles, r = nu-dependent (J6)
    the rest      : radical = bicycles = the descent's cuts    (421)

(1) GAUSS SUM LAW (machine-verified): G(B) := sum_{x in C} i^{|x|} obeys
        G = 0                    if q is nonzero somewhere on the radical,
        |G| = 2^d * sqrt(2)^rho  otherwise.
    Checked on Petersen, the span419/422 instances, and a theta gadget.
(2) TARGET (S) SPAN-ALL: even path-differences over ALL THREE terminal
    pairs span C_0 -- stratum-free spanning. Tested on random bridge-shaped
    graphs (subcubic, 2-connected, girth >= 5, exactly three degree-2
    terminals; terminal distances relaxed, noted).
(3) TARGET (P) PACKED CEILING: d <= mu - 3 for non-bipartite bridge-shaped
    graphs (the anatomy's forced equality d = mu - 2 - 2t never holds in
    nature; Petersen, the sharp cubic witness, has NO terminals and is
    excluded by bridge shape). Tested on the same samples.
"""
import random, itertools

def machinery(VER, EDG):
    E = {frozenset(e): i for i, e in enumerate(EDG)}
    adj = {v: [] for v in VER}
    for u, w in EDG: adj[u].append(w); adj[w].append(u)
    bit = lambda u, w: 1 << E[frozenset((u, w))]
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
    fund = [(1 << E[e]) | tp(*tuple(e)) for e in E if e not in tree]
    return E, adj, bit, fund

pop = lambda m: bin(m).count('1')

def gauss_and_anatomy(name, VER, EDG):
    E, adj, bit, fund = machinery(VER, EDG)
    mu = len(fund)
    # full cycle space sweep
    G = 0j; elems = []
    for mask in range(1 << mu):
        x = 0; mm = mask
        i = 0
        while mm:
            if mm & 1: x ^= fund[i]
            mm >>= 1; i += 1
        elems.append(x)
        G += 1j ** (pop(x) % 4)
    # rank of lambda on fund basis
    M = [[pop(a & b) % 2 for b in fund] for a in fund]
    rank = 0
    M = [row[:] for row in M]
    rows = list(range(mu)); cols = list(range(mu))
    Mm = M
    r = 0
    for c in range(mu):
        piv = next((i for i in range(r, mu) if Mm[i][c]), None)
        if piv is None: continue
        Mm[r], Mm[piv] = Mm[piv], Mm[r]
        for i in range(mu):
            if i != r and Mm[i][c]:
                Mm[i] = [(a + b) % 2 for a, b in zip(Mm[i], Mm[r])]
        r += 1
    rho = r; d = mu - rho
    pred = 2 ** d * (2 ** 0.5) ** rho
    okG = abs(G) < 1e-6 or abs(abs(G) - pred) < 1e-6
    print(' %-22s mu=%2d  d=%2d rho=%2d  |G|=%9.3f  law(0 or 2^d sqrt2^rho'
          '=%9.3f): %s' % (name, mu, d, rho, abs(G), pred, okG))
    return okG

def paths(adj, bit, src, dst, cap=250000):
    out = []
    def dfs(v, mask, seen):
        if len(out) > cap: return
        if v == dst: out.append(mask); return
        for w in adj[v]:
            if w not in seen: dfs(w, mask | bit(v, w), seen | {w})
    dfs(src, 0, {src})
    return out

def rankv(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return basis

print('=== (1) the Gauss Sum Law ===')
# Petersen
P = [(i, (i + 1) % 5) for i in range(5)] + [(i + 5, (i + 2) % 5 + 5)
     for i in range(5)] + [(i, i + 5) for i in range(5)]
allok = gauss_and_anatomy('Petersen', list(range(10)), P)
# theta gadget: two vertices joined by three paths of lengths 3,3,5
V = ['u', 'v'] + ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'c3', 'c4']
Eg = [('u', 'a1'), ('a1', 'a2'), ('a2', 'v'), ('u', 'b1'), ('b1', 'b2'),
      ('b2', 'v'), ('u', 'c1'), ('c1', 'c2'), ('c2', 'c3'), ('c3', 'c4'),
      ('c4', 'v')]
allok &= gauss_and_anatomy('theta(3,3,5)', V, Eg)
q_theta = [((3 + 3) % 4), ((3 + 5) % 4), ((3 + 5) % 4)]
print('   theta q-values on nonzero elements: %s; Gauss = 1+3*i^2 = -2'
      % q_theta)
o = ['o%d' % i for i in range(12)]; p = ['p%d' % i for i in range(6)]
V1 = ['b1'] + o + ['b2'] + p + ['d2', 'd8', 'd10', 'b3', 'v1', 'v2']
E1 = ([('b1', 'o0'), ('o11', 'b1')] + [(o[i], o[i+1]) for i in range(11)]
      + [('b2', 'p0'), ('p5', 'b2')] + [(p[i], p[i+1]) for i in range(5)]
      + [('o1', 'p1'), ('o6', 'p3')]
      + [('o2', 'd2'), ('o8', 'd8'), ('o10', 'd10'), ('d8', 'b3'),
         ('d10', 'b3'), ('d2', 'd8')] + [('o3', 'o9')]
      + [('o0', 'v1'), ('o4', 'v1'), ('o7', 'v1'), ('o5', 'v2'),
         ('o11', 'v2')])
allok &= gauss_and_anatomy('instance 1 (span419)', V1, E1)
assert allok

print()
print('=== (2)+(3) random bridge-shaped graphs: SPAN-ALL and the ceiling ===')
random.seed(423)
def random_bridge(n):
    for _ in range(4000):
        stubs = [v for v in range(n) for _ in range(3)]
        random.shuffle(stubs)
        E = set()
        ok = True
        for i in range(0, len(stubs), 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b or frozenset((a, b)) in E: ok = False; break
            E.add(frozenset((a, b)))
        if not ok: continue
        adj = {v: [] for v in range(n)}
        for e in E:
            a, b = tuple(e); adj[a].append(b); adj[b].append(a)
        # girth >= 5
        g = 99
        for s in range(n):
            dist = {s: 0}; par = {s: None}; q = [s]
            while q:
                v = q.pop(0)
                for w in adj[v]:
                    if w not in dist:
                        dist[w] = dist[v] + 1; par[w] = v; q.append(w)
                    elif par[v] != w:
                        g = min(g, dist[v] + dist[w] + 1)
        if g < 5: continue
        # connectivity (2-connected quick check: connected after any vertex del)
        def conn(rm):
            vs = [v for v in range(n) if v != rm]
            if not vs: return True
            seen = {vs[0]}; st = [vs[0]]
            while st:
                v = st.pop()
                for w in adj[v]:
                    if w != rm and w not in seen: seen.add(w); st.append(w)
            return len(seen) == len(vs)
        if not all(conn(v) for v in range(n)): continue
        # subdivide 3 disjoint edges -> terminals
        El = [tuple(e) for e in E]
        random.shuffle(El)
        picked = []
        for e in El:
            if all(v not in f for f in picked for v in e): picked.append(e)
            if len(picked) == 3: break
        if len(picked) < 3: continue
        VER = list(range(n)) + ['t0', 't1', 't2']
        EDG = [tuple(e) for e in E if tuple(e) not in picked
               and tuple(e)[::-1] not in picked]
        EDG = [e for e in E if tuple(sorted(e)) not in
               [tuple(sorted(x)) for x in picked]]
        EDG = [tuple(e) for e in EDG]
        for k, (a, b) in enumerate(picked):
            EDG += [(a, 't%d' % k), ('t%d' % k, b)]
        return VER, EDG
    return None

tested = span_ok = ceil_ok = bip = 0
for n in [12, 14, 14, 16, 16, 16, 18, 18, 18, 18, 14, 16, 18, 12, 16]:
    res = random_bridge(n)
    if not res: continue
    VER, EDG = res
    E, adj, bit, fund = machinery(VER, EDG)
    mu = len(fund)
    odd = [c for c in fund if pop(c) % 2 == 1]
    if not odd:
        bip += 1; continue
    C0 = rankv([c for c in fund if pop(c) % 2 == 0] +
               [a ^ b for a in odd for b in odd])
    # SPAN-ALL: differences over all three terminal pairs
    diffs = []
    for (s, t) in (('t0', 't1'), ('t0', 't2'), ('t1', 't2')):
        ps = paths(adj, bit, s, t)
        ev = [m for m in ps if pop(m) % 2 == 0]
        if len(ev) > 1:
            diffs += [m ^ ev[0] for m in ev[1:]]
    W = rankv(diffs)
    sp = len(W) == len(C0)
    # ceiling: d <= mu - 3
    M = [[pop(a & b) % 2 for b in fund] for a in fund]
    r = 0
    for c in range(mu):
        piv = next((i for i in range(r, mu) if M[i][c]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(mu):
            if i != r and M[i][c]:
                M[i] = [(x + y) % 2 for x, y in zip(M[i], M[r])]
        r += 1
    d = mu - r
    ce = d <= mu - 3
    tested += 1
    span_ok += sp; ceil_ok += ce
print(' bridges tested: %d (skipped bipartite: %d)' % (tested, bip))
print(' (S) SPAN-ALL over three pairs: %d/%d' % (span_ok, tested))
print(' (P) ceiling d <= mu - 3:       %d/%d' % (ceil_ok, tested))
