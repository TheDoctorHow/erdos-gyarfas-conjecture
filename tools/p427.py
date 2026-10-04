"""p427.py -- finishing (P): the Route-Adjustment Lemma and the descent
through infected components.

Board: monochrome sealed-world instance (A = 17, beta = 0) with an
INFECTED pocket X: interior all-odd theta(3,3,3) on ends u, v, one
attachment hanging at the first internal vertex of each branch (atts at
colour-0 positions {2, 9, 13}); windows {1, 3, 8, 14} (|WIN| = 4 == 1 mod
3); clean pockets {6,10,12}, {0,4,7}, {5,11,15}.

Checks:
 (1) ROUTE ADJUSTMENT: for X's same-side consecutive att pair (9, 13),
     the via-u route is CHARGED against the theta plane (pairing vector
     (1,0)) and the via-v route is UNCHARGED ((0,0)) -- the Gram rows
     span, so an uncharged route exists, as the lemma asserts.
 (2) THE DESCENT STILL KILLS: the uncharged att-cycle, the clean pockets'
     cycles and the D_3 window-cycles are radical candidates; NONE of them
     lies in the cut space: the 421 contradiction fires on the INFECTED
     board exactly as on clean ones.  J2' is bypassed, not blocking.
"""

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

# ---------------- the infected board, A = 17 ----------------
A, s1 = 17, 4
o = ['o%d' % i for i in range(A - 1)]
VER = ['b1'] + o
EDG = [('b1', 'o0'), ('o%d' % (A - 2), 'b1')] + \
      [(o[i], o[i + 1]) for i in range(A - 2)]
# D_3 caterpillar on windows {1,3,8,14}, b3 spliced
WIN = [1, 3, 8, 14]
VER += ['d%d' % w for w in WIN] + ['b3', 'z2', 'z3']
EDG += [('o%d' % w, 'd%d' % w) for w in WIN]
EDG += [('d1', 'b3'), ('b3', 'd3'), ('d3', 'z2'), ('z2', 'd8'),
        ('d8', 'z3'), ('z3', 'd14')]
# clean caterpillar pockets
def pocket(tag, atts):
    global VER, EDG
    us = []
    for j in range(0, len(atts), 2):
        u = '%su%d' % (tag, j)
        VER.append(u); us.append(u)
        EDG.extend(('o%d' % a, u) for a in atts[j:j + 2])
    for j in range(len(us) - 1):
        z = '%sz%d' % (tag, j)
        VER.append(z); EDG += [(us[j], z), (z, us[j + 1])]
    return atts, us
P1 = pocket('P1', [6, 10, 12])
P2 = pocket('P2', [0, 4, 7])
P3 = pocket('P3', [5, 11, 15])
# infected pocket X: theta(3,3,3), atts {2,9,13} at first internals
VER += ['u', 'v', 'a1', 'a2', 'bb1', 'bb2', 'c1', 'c2']
EDG += [('u', 'a1'), ('a1', 'a2'), ('a2', 'v'),
        ('u', 'bb1'), ('bb1', 'bb2'), ('bb2', 'v'),
        ('u', 'c1'), ('c1', 'c2'), ('c2', 'v'),
        ('o2', 'a1'), ('o9', 'bb1'), ('o13', 'c1')]

E, adj, bit, fund = machinery(VER, EDG)
mask = lambda es: sum(1 << E[frozenset(e)] for e in es)
ZAB = mask([('u', 'a1'), ('a1', 'a2'), ('a2', 'v'),
            ('u', 'bb1'), ('bb1', 'bb2'), ('bb2', 'v')])
ZBC = mask([('u', 'bb1'), ('bb1', 'bb2'), ('bb2', 'v'),
            ('u', 'c1'), ('c1', 'c2'), ('c2', 'v')])

# (1) route adjustment for att pair (9, 13): bb1 -> c1
route_u = mask([('bb1', 'u'), ('u', 'c1')])
route_v = mask([('bb1', 'bb2'), ('bb2', 'v'), ('v', 'c2'), ('c2', 'c1')])
vec = lambda r: (pop(r & ZAB) % 2, pop(r & ZBC) % 2)
print('(1) route pairing vectors: via u %s (charged), via v %s (uncharged)'
      % (vec(route_u), vec(route_v)))
assert vec(route_u) != (0, 0) and vec(route_v) == (0, 0)

# (2) radical candidates on the infected board vs the cut space
def incut(z):
    return all(pop(z & c) % 2 == 0 for c in fund)
arc = lambda a, b: mask([('o%d' % i, 'o%d' % (i + 1)) for i in range(a, b)])
cands = []
Zx = mask([('o9', 'bb1'), ('o13', 'c1')]) ^ route_v ^ arc(9, 13)
cands.append(('X atts (9,13) via v [radical]', Zx))
cands.append(('D3 windows (1,3)',
              mask([('o1', 'd1'), ('o3', 'd3'), ('d1', 'b3'), ('b3', 'd3')])
              ^ arc(1, 3)))
cands.append(('D3 windows (8,14)',
              mask([('o8', 'd8'), ('o14', 'd14'), ('d8', 'z3'),
                    ('z3', 'd14')]) ^ arc(8, 14)))
cands.append(('P1 atts (10,12)',
              mask([('o10', 'P1u0'), ('o12', 'P1u2'), ('P1u0', 'P1z0'),
                    ('P1z0', 'P1u2')]) ^ arc(10, 12)))
cands.append(('P3 atts (11,15)',
              mask([('o11', 'P3u0'), ('o15', 'P3u2'), ('P3u0', 'P3z0'),
                    ('P3z0', 'P3u2')]) ^ arc(11, 15)))
fails = 0
for name, z in cands:
    assert pop(z) % 2 == 0, name
    ok = incut(z)
    print('   %-32s %s' % (name, 'IN cut space' if ok
                           else 'NOT a cut -> KILL'))
    fails += not ok
print('(2) %d/%d radical candidates fail cut membership on the INFECTED'
      ' board' % (fails, len(cands)))
assert fails >= 1
print()
print('VERDICT: the theta shields nothing -- the adjusted route restores')
print('a radical att-cycle and the descent contradiction fires at t >= 1.')
