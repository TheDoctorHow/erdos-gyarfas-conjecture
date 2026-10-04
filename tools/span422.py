"""span422.py -- closing the joints J1-J3 sweepingly (audit for 422).

(1) {1,2}-SIGNATURE LAW: every b_1-b_2 path satisfies
        |P| == chi(d_1) + chi(d_2) + [f in P]   (mod 2),
    d_1, d_2 the doors. Exhaustive over three instances (b_3-interior
    passages included). The parity repair is a DOOR SWAP at whichever end
    is away from the target edge: no mid-path surgery, no collisions.
(2) J1 CLOSED: W_12 := span{P xor P' : P, P' even b_1-b_2 paths} = C_0,
    machine-verified end to end on the two span419 canonical instances.
(3) SWAP QUADRUPLE: the algebraic closure P_{u',y} = P_{u',x} xor P_{u,y}
    xor P_{u,x} -- all four signature paths present over a common course,
    so 222.2's fixed-base lemma applies with zero routing hypotheses.
(4) THE WEAVE: instance 3 carries two cross-colour INTERLEAVED pockets
    (X at {7,11,13} colour 0, X' at {6,10,12} colour 1; pairs (7,11) and
    (10,12) interleave with odd overlap [10,11]). The four paths P_0,
    P_0+Z, P_0+Z', P_0+Z+Z' are all genuine even paths, so under (SC)
    222.2 forces |Z cap Z'| even -- but it is 1. Boards with cross-colour
    interleaving are EXECUTED.
"""

def machinery(VER, EDG):
    E = {}
    for e in EDG: E[frozenset(e)] = len(E)
    adj = {v: [] for v in VER}
    for u, w in EDG: adj[u].append(w); adj[w].append(u)
    bit = lambda u, w: 1 << E[frozenset((u, w))]
    def paths(src, dst):
        out = []
        def dfs(v, mask, seen):
            if v == dst: out.append(mask); return
            for w in adj[v]:
                if w not in seen: dfs(w, mask | bit(v, w), seen | {w})
        dfs(src, 0, {src}); return out
    tree = set(); seen = {VER[0]}; st = [VER[0]]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen: seen.add(w); tree.add(frozenset((v, w))); st.append(w)
    par = {VER[0]: None}; st = [VER[0]]
    while st:
        v = st.pop()
        for w in adj[v]:
            if frozenset((v, w)) in tree and w not in par: par[w] = v; st.append(w)
    def tp(a, b):
        m = 0; anc = set(); x = a
        while x is not None: anc.add(x); x = par[x]
        y = b
        while y not in anc: m |= bit(y, par[y]); y = par[y]
        x = a
        while x != y: m |= bit(x, par[x]); x = par[x]
        return m
    fund = [(1 << E[e]) | tp(*tuple(e)) for e in E if e not in tree]
    return E, adj, bit, paths, fund

pop = lambda m: bin(m).count('1')
def rank(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return basis
def inspan(v, basis):
    for b in sorted(basis, reverse=True): v = min(v, v ^ b)
    return v == 0

def colour(VER, EDG, b1, b2, f):
    adj = {v: [] for v in VER}
    for u, w in EDG: adj[u].append(w); adj[w].append(u)
    Hv = [v for v in VER if v not in (b1, b2)]
    chi = {Hv[0]: 0}; st = [Hv[0]]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w in (b1, b2) or frozenset((v, w)) == frozenset(f): continue
            if w not in chi: chi[w] = chi[v] ^ 1; st.append(w)
            else: assert chi[w] != chi[v]
    return chi

def run(name, VER, EDG, b1, b2, f, weave=None):
    E, adj, bit, paths, fund = machinery(VER, EDG)
    chi = colour(VER, EDG, b1, b2, f)
    fbit = 1 << E[frozenset(f)]
    P = paths(b1, b2)
    ok = 0
    for m in P:
        d1 = [w for w in adj[b1] if m & bit(b1, w)][0]
        d2 = [w for w in adj[b2] if m & bit(b2, w)][0]
        pred = (chi[d1] + chi[d2] + (1 if m & fbit else 0)) % 2
        if pop(m) % 2 == pred: ok += 1
    print(' %s: signature law %d/%d paths conform' % (name, ok, len(P)))
    assert ok == len(P)
    ev = [m for m in P if pop(m) % 2 == 0]
    C0 = rank([c for c in fund if pop(c) % 2 == 0] +
              [a ^ b for a in fund for b in fund
               if pop(a) % 2 == 1 and pop(b) % 2 == 1])
    W = rank([m ^ ev[0] for m in ev[1:]])
    good = all(inspan(w, C0) for w in W) and len(W) == len(C0)
    print('   W_12: dim %d = dim C_0 %d, equal: %s   (mu = %d)'
          % (len(W), len(C0), good, len(fund)))
    assert good
    # swap quadruple: four paths with P4 = P1^P2^P3, doors u/u' x x/y
    u1, u2 = adj[b1]; x1, x2 = adj[b2]
    evset = set(ev)
    sig = {}
    for m in ev:
        d1 = u1 if m & bit(b1, u1) else u2
        d2 = x1 if m & bit(b2, x1) else x2
        sig.setdefault((d1, d2), []).append(m)
    found = False
    for p1 in sig.get((u1, x1), [])[:80]:
        for p2 in sig.get((u2, x1), [])[:80]:
            for p3 in sig.get((u1, x2), [])[:80]:
                if p1 ^ p2 ^ p3 in evset:
                    found = True; break
            if found: break
        if found: break
    print('   swap quadruple (u,x),(u\',x),(u,y),(u\',y) with closure: %s'
          % found)
    assert found
    if weave:
        Z1, Z2 = weave
        m1 = 0
        for e in Z1: m1 |= 1 << E[frozenset(e)]
        m2 = 0
        for e in Z2: m2 |= 1 << E[frozenset(e)]
        quad = [q for q in (lambda base: [(b0, b0 ^ m1, b0 ^ m2,
                b0 ^ m1 ^ m2) for b0 in base])(ev)
                if all(x in evset for x in q)]
        print('   WEAVE: |Z cap Z\'| = %d (odd: %s); full four-path'
              ' quadruples found: %d  -> (SC) EXECUTES the board'
              % (pop(m1 & m2), pop(m1 & m2) % 2 == 1, len(quad)))
        assert pop(m1 & m2) % 2 == 1 and quad

# instance 1 (span419): A = 13 canonical
o = ['o%d' % i for i in range(12)]; p = ['p%d' % i for i in range(6)]
V1 = ['b1'] + o + ['b2'] + p + ['d2', 'd8', 'd10', 'b3', 'v1', 'v2']
E1 = ([('b1', 'o0'), ('o11', 'b1')] + [(o[i], o[i+1]) for i in range(11)]
      + [('b2', 'p0'), ('p5', 'b2')] + [(p[i], p[i+1]) for i in range(5)]
      + [('o1', 'p1'), ('o6', 'p3')]
      + [('o2', 'd2'), ('o8', 'd8'), ('o10', 'd10'),
         ('d8', 'b3'), ('d10', 'b3'), ('d2', 'd8')]
      + [('o3', 'o9')]
      + [('o0', 'v1'), ('o4', 'v1'), ('o7', 'v1'), ('o5', 'v2'),
         ('o11', 'v2')])
run('instance 1 (A=13)', V1, E1, 'b1', 'b2', ('o4', 'o5'))

# instance 2 (span419): A = 11 canonical
o = ['o%d' % i for i in range(10)]; p = ['p%d' % i for i in range(4)]
V2 = ['b1'] + o + ['b2'] + p + ['d3', 'd8', 'b3', 'v1', 'v2']
E2 = ([('b1', 'o0'), ('o9', 'b1')] + [(o[i], o[i+1]) for i in range(9)]
      + [('b2', 'p0'), ('p3', 'b2')] + [(p[i], p[i+1]) for i in range(3)]
      + [('o1', 'p1'), ('o5', 'p3')]
      + [('o3', 'd3'), ('o8', 'd8'), ('d3', 'b3'), ('d8', 'b3')]
      + [('o0', 'v1'), ('o4', 'v1'), ('o7', 'v1'), ('o2', 'v2'),
         ('o6', 'v2'), ('o9', 'v2')])
run('instance 2 (A=11)', V2, E2, 'b1', 'b2', ('o6', 'o7'))

# instance 3: A = 15 with cross-colour INTERLEAVED pockets (the weave board)
o = ['o%d' % i for i in range(14)]; p = ['p%d' % i for i in range(6)]
V3 = (['b1'] + o + ['b2'] + p
      + ['d1', 'd3', 'd8', 'b3', 'm']
      + ['xu1', 'xz', 'xu2', 'yu1', 'yz', 'yu2', 'pu1', 'pz', 'pu2',
         'c1', 'c2'])
E3 = ([('b1', 'o0'), ('o13', 'b1')] + [(o[i], o[i+1]) for i in range(13)]
      + [('b2', 'p0'), ('p5', 'b2')] + [(p[i], p[i+1]) for i in range(5)]
      + [('o2', 'c1'), ('c1', 'p1'), ('o5', 'c2'), ('c2', 'p3')]
      + [('o1', 'd1'), ('o3', 'd3'), ('o8', 'd8'),
         ('d1', 'b3'), ('b3', 'd3'), ('d3', 'm'), ('m', 'd8')]
      + [('o7', 'xu1'), ('o11', 'xu1'), ('xu1', 'xz'), ('xz', 'xu2'),
         ('o13', 'xu2')]
      + [('o6', 'yu1'), ('o10', 'yu1'), ('yu1', 'yz'), ('yz', 'yu2'),
         ('o12', 'yu2')]
      + [('o0', 'pu1'), ('o4', 'pu1'), ('pu1', 'pz'), ('pz', 'pu2'),
         ('o9', 'pu2')])
ZX = [('o7', 'xu1'), ('o11', 'xu1')] + [('o%d' % i, 'o%d' % (i+1))
                                        for i in range(7, 11)]
ZY = [('o10', 'yu1'), ('yu1', 'yz'), ('yz', 'yu2'), ('yu2', 'o12')] \
     + [('o%d' % i, 'o%d' % (i+1)) for i in range(10, 12)]
run('instance 3 (A=15, weave)', V3, E3, 'b1', 'b2', ('o4', 'o5'),
    weave=(ZX, ZY))
print()
print('J1: W_12 = C_0 with end-local repairs -- no arc-ordering cases.')
print('J2: swap closure + parity classification + weave execution.')
print('J3: monotonicity (see 422.4) -- audited by grep, no machine part.')
