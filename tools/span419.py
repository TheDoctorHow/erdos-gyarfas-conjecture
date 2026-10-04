"""span419.py -- target (i): COVERAGE ==> SPANNING in the sealed canonical world.

Two synthetic instances of the 400.3 canonical form (O_1 an odd cycle through
b_1 carrying the corridor edge f; O_2 a disjoint odd cycle through b_2; H =
B - {b_1,b_2} - f bipartite; D_3 with b_3 straddling; two transitions, forced
by 2-connectivity; chords/pockets covering the remaining O_1 positions).

Checks per instance:
 1. geometry: subcubic, H bipartite, colour laws 400.1 (i)-(iii), b_3 opposite.
 2. DOOR PARITY LAWS (exhaustive over all simple connecting paths):
    {1,3}: parity == 1 + chi(door) + chi(b_3) + [uses f] + [passes b_2]
      in particular pure-H {1,3}-paths ALL even; f-only ODD (G1 recovered and
      extended); b_2-only ODD; f-and-b_2 EVEN (the two-transition class G1
      never saw).
    {2,3}: parity == 1 + chi(door) + chi(b_3) + [uses f].
 3. THE SPANNING THEOREM: W := span{ P xor P_0 : P even b_1-b_3 path } has
    dim W = mu - 1 and W = C_0 (the even-weight cycle subspace).
 4. THE DECOMPOSITION: C_0 = C(H) (+) <gamma_1> (+) <gamma_2 + gamma_3>:
    on C_0 the f-coordinate equals the b_2-coordinate, and both b-coordinates
    are realized by explicit differences.
"""
import itertools

def solve(name, VER, EDG, b1, b2, b3, f_edge, b2_edges, doors13, doors23):
    print('--- instance %s ---' % name)
    V = {v: i for i, v in enumerate(VER)}
    E = {}
    for e in EDG:
        E[frozenset(e)] = len(E)
    nE = len(E)
    adj = {v: [] for v in VER}
    for e in EDG:
        u, w = e
        adj[u].append(w); adj[w].append(u)
    assert all(len(adj[v]) <= 3 for v in VER), 'subcubic'
    assert all(len(adj[v]) == 2 for v in (b1, b2, b3)), 'terminals degree 2'

    # H = B - b1 - b2 - f : bipartite 2-colouring
    Hverts = [v for v in VER if v not in (b1, b2)]
    Hadj = {v: [w for w in adj[v] if w not in (b1, b2)
               and frozenset((v, w)) != frozenset(f_edge)] for v in Hverts}
    chi = {Hverts[0]: 0}
    stack = [Hverts[0]]
    while stack:
        v = stack.pop()
        for w in Hadj[v]:
            if w not in chi:
                chi[w] = chi[v] ^ 1; stack.append(w)
            else:
                assert chi[w] != chi[v], 'H not bipartite'
    assert len(chi) == len(Hverts), 'H connected'
    g, h = f_edge
    u1, u2 = [w for w in adj[b1]]
    x1, x2 = [w for w in adj[b2]]
    print(' colour laws: chi(g)=chi(h): %s; chi(u)=chi(u\'): %s; '
          'chi(x)!=chi(y): %s; b3 opposite quadruple: %s'
          % (chi[g] == chi[h], chi[u1] == chi[u2], chi[x1] != chi[x2],
             chi[b3] != chi[u1]))
    assert chi[g] == chi[h] and chi[u1] == chi[u2] and chi[x1] != chi[x2]

    def bit(u, w):
        return 1 << E[frozenset((u, w))]
    fbit = 1 << E[frozenset(f_edge)]
    b2bits = 0
    for e in b2_edges:
        b2bits |= 1 << E[frozenset(e)]

    def paths(src, dst, forbid=()):
        out = []
        def dfs(v, mask, seen):
            if v == dst:
                out.append(mask); return
            for w in adj[v]:
                if w in seen or w in forbid: continue
                dfs(w, mask | bit(v, w), seen | {w})
        dfs(src, 0, {src})
        return out

    def plen(mask): return bin(mask).count('1')

    # 2. door parity laws
    for (src, dst, doors, tag) in ((b1, b3, doors13, '{1,3}'),
                                   (b2, b3, doors23, '{2,3}')):
        ok = tot = 0
        for m in paths(src, dst):
            door = next(d for d in doors if (m >> E[frozenset((src, d))]) & 1
                        if True) if False else None
            door = [d for d in doors if m & bit(src, d)][0]
            usef = 1 if m & fbit else 0
            useb2 = 1 if (src != b2 and (m & b2bits) == b2bits) else 0
            pred = (1 + chi[door] + chi[b3] + usef + useb2) % 2
            tot += 1
            if plen(m) % 2 == pred: ok += 1
        print(' door parity law %s: %d/%d paths conform' % (tag, ok, tot))
        assert ok == tot

    # linear algebra over F_2
    def rank(vectors):
        basis = []
        for v in vectors:
            for b in basis:
                v = min(v, v ^ b)
            if v: basis.append(v)
        return basis
    def in_span(v, basis):
        for b in sorted(basis, reverse=True):
            v = min(v, v ^ b)
        return v == 0

    # cycle space C and even part C_0
    tree = set(); seen = {VER[0]}
    stack = [VER[0]]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w); tree.add(frozenset((v, w))); stack.append(w)
    def treepath(a, bb):
        par = {a: None}
        st = [a]
        while st:
            v = st.pop()
            for w in adj[v]:
                if frozenset((v, w)) in tree and w not in par:
                    par[w] = v; st.append(w)
        m = 0; v = bb
        while par[v] is not None:
            m |= bit(v, par[v]); v = par[v]
        return m
    fund = []
    for e in E:
        if e not in tree:
            u, w = tuple(e)
            fund.append((1 << E[e]) | treepath(u, w))
    mu = len(fund)
    C = rank(fund)
    C0 = rank([c for c in fund if plen(c) % 2 == 0] +
              [c1 ^ c2 for c1 in fund for c2 in fund
               if plen(c1) % 2 == 1 and plen(c2) % 2 == 1])
    print(' mu = %d; dim C = %d; dim C_0 = %d (= mu - 1: %s)'
          % (mu, len(C), len(C0), len(C0) == mu - 1))

    # 3. the spanning check
    evens = [m for m in paths(b1, b3) if plen(m) % 2 == 0]
    P0 = evens[0]
    W = rank([m ^ P0 for m in evens[1:]])
    inC0 = all(in_span(w, C0) for w in W)
    print(' SPANNING: dim W = %d; W subseteq C_0: %s; W = C_0: %s'
          % (len(W), inC0, inC0 and len(W) == len(C0)))
    assert inC0 and len(W) == len(C0)

    # 4. decomposition: on C_0 the f-coordinate == b_2-pair-coordinate;
    #    realizers for the b_1- and (b_2,f)-directions exist among differences
    b1bits = bit(b1, u1) | bit(b1, u2)
    for c in C0:
        cf = 1 if c & fbit else 0
        cb2 = 1 if (c & b2bits) == b2bits else 0
        assert (c & b2bits) in (0, b2bits) and cf == cb2
    d_b1 = any(((m ^ P0) & b1bits) == b1bits for m in evens)
    d_b2f = any((m ^ P0) & fbit and ((m ^ P0) & b2bits) == b2bits
                for m in evens)
    print(' decomposition: f-coord == b_2-coord on all of C_0: True;'
          ' b_1-direction realized: %s; (b_2+f)-direction realized: %s'
          % (d_b1, d_b2f))
    assert d_b1 and d_b2f
    print(' instance %s: ALL CHECKS PASS' % name)
    print()

# ---------- instance 1: A = 13, s_1 = 4, windows {2,8,10}, chord, 2 pockets
o = ['o%d' % i for i in range(12)]
p = ['p%d' % i for i in range(6)]
VER1 = ['b1'] + o + ['b2'] + p + ['d2', 'd8', 'd10', 'b3', 'v1', 'v2']
EDG1 = ([('b1', 'o0'), ('o11', 'b1')]
        + [(o[i], o[i + 1]) for i in range(11)]              # incl. f = o4-o5
        + [('b2', 'p0'), ('p5', 'b2')]
        + [(p[i], p[i + 1]) for i in range(5)]
        + [('o1', 'p1'), ('o6', 'p3')]                       # two transitions
        + [('o2', 'd2'), ('o8', 'd8'), ('o10', 'd10'),
           ('d8', 'b3'), ('d10', 'b3'), ('d2', 'd8')]        # D_3
        + [('o3', 'o9')]                                     # chord
        + [('o0', 'v1'), ('o4', 'v1'), ('o7', 'v1'),
           ('o5', 'v2'), ('o11', 'v2')])                     # pockets
solve('1 (A=13)', VER1, EDG1, 'b1', 'b2', 'b3', ('o4', 'o5'),
      [('b2', 'p0'), ('p5', 'b2')], ('o0', 'o11'), ('p0', 'p5'))

# ---------- instance 2: A = 11, s_1 = 6, windows {3,8}, 2 pockets
o = ['o%d' % i for i in range(10)]
p = ['p%d' % i for i in range(4)]
VER2 = ['b1'] + o + ['b2'] + p + ['d3', 'd8', 'b3', 'v1', 'v2']
EDG2 = ([('b1', 'o0'), ('o9', 'b1')]
        + [(o[i], o[i + 1]) for i in range(9)]               # incl. f = o6-o7
        + [('b2', 'p0'), ('p3', 'b2')]
        + [(p[i], p[i + 1]) for i in range(3)]
        + [('o1', 'p1'), ('o5', 'p3')]                       # two transitions
        + [('o3', 'd3'), ('o8', 'd8'), ('d3', 'b3'), ('d8', 'b3')]
        + [('o0', 'v1'), ('o4', 'v1'), ('o7', 'v1'),
           ('o2', 'v2'), ('o6', 'v2'), ('o9', 'v2')])        # pockets
solve('2 (A=11)', VER2, EDG2, 'b1', 'b2', 'b3', ('o6', 'o7'),
      [('b2', 'p0'), ('p3', 'b2')], ('o0', 'o9'), ('p0', 'p3'))

print('TARGET (i) INSTANCE VERDICT: spanning holds on every instance tested.')
