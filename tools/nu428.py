"""nu428.py -- the portable demolition: the engine on a nu = 1 world with
no corridor, no template, no seal.

Board: a single odd cycle Lambda (length 11) carrying terminal b_1 at
position 0; every other position carries its third edge (coverage = cubic
ambient); three middle components, b_3 inside X1, b_2 inside X2, all
bipartite; NO second odd cycle, NO corridor edge, NO windows/feet.

Checks:
 (1) MONO-PAIR PIGEONHOLE: every >= 3-attachment component has an
     attachment pair at even Lambda-distance (2-colouring pigeonhole).
 (2) the even-pair att-cycles (even arc => orthogonal to Lambda, i.e.
     radical candidates under the nu = 1 anatomy) are NOT cuts: the
     descent contradiction fires with r = 1, no opposed pair, no seal.
 (3) odd-pair att-cycles are Lambda-CHARGED (the r = 1 escape hatch) --
     present on the board, and irrelevant: the descent never needs them.
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

N = 11
L = ['L%d' % i for i in range(N)]          # L0 = b_1 (degree 2, on Lambda)
VER = L[:]
EDG = [(L[i], L[(i + 1) % N]) for i in range(N)]
# X1: atts {1,4,7,9}, contains b_3
VER += ['w1', 'w4', 'w7', 'w9', 'b3', 't1', 't2']
EDG += [('L1', 'w1'), ('L4', 'w4'), ('L7', 'w7'), ('L9', 'w9'),
        ('w1', 'b3'), ('b3', 't1'), ('t1', 'w4'), ('w4', 'w7'),
        ('w7', 't2'), ('t2', 'w9')]
# X2: atts {2,5,8}, contains b_2
VER += ['w2', 'w5', 'w8', 'b2', 'u']
EDG += [('L2', 'w2'), ('L5', 'w5'), ('L8', 'w8'),
        ('w2', 'b2'), ('b2', 'u'), ('u', 'w5'), ('w5', 'w8')]
# X3: atts {3,6,10}
VER += ['w3', 'w6', 'w10', 'vv']
EDG += [('L3', 'w3'), ('L6', 'w6'), ('L10', 'w10'),
        ('w3', 'w6'), ('w6', 'vv'), ('vv', 'w10')]

E, adj, bit, fund = machinery(VER, EDG)
mask = lambda es: sum(1 << E[frozenset(e)] for e in es)
Lam = mask([(L[i], L[(i + 1) % N]) for i in range(N)])
arc = lambda a, b: mask([(L[i], L[i + 1]) for i in range(a, b)])
incut = lambda z: all(pop(z & c) % 2 == 0 for c in fund)

# (1) pigeonhole
comps = {'X1': [1, 4, 7, 9], 'X2': [2, 5, 8], 'X3': [3, 6, 10]}
for name, atts in comps.items():
    evenpairs = [(a, b) for i, a in enumerate(atts) for b in atts[i + 1:]
                 if (b - a) % 2 == 0]
    print('(1) %s atts %s: even-distance pairs %s' % (name, atts, evenpairs))
    assert evenpairs

# (2) radical att-cycles (even arcs) vs the cut space
paths = {('X1', 7, 9): [('w7', 't2'), ('t2', 'w9')],
         ('X1', 1, 7): [('w1', 'b3'), ('b3', 't1'), ('t1', 'w4'),
                        ('w4', 'w7')],
         ('X2', 2, 8): [('w2', 'b2'), ('b2', 'u'), ('u', 'w5'),
                        ('w5', 'w8')],
         ('X3', 6, 10): [('w6', 'vv'), ('vv', 'w10')]}
kills = 0
for (name, a, b), interior in paths.items():
    z = mask([('L%d' % a, 'w%d' % a), ('L%d' % b, 'w%d' % b)]) \
        ^ mask(interior) ^ arc(a, b)
    assert pop(z) % 2 == 0 and pop(z & Lam) % 2 == 0, 'not radical-shaped'
    ok = incut(z)
    print('(2) %s att-cycle (%d,%d): %s' % (name, a, b,
          'IN cut space' if ok else 'NOT a cut -> KILL'))
    kills += not ok
assert kills >= 1
# (3) a charged odd-pair cycle, for the record
z58 = mask([('L5', 'w5'), ('L8', 'w8'), ('w5', 'w8')]) ^ arc(5, 8)
print('(3) X2 odd pair (5,8): lambda(., Lambda) = %d  -> CHARGED, unused'
      % (pop(z58 & Lam) % 2))
assert pop(z58 & Lam) % 2 == 1
print()
print('%d/%d radical candidates fail cut membership on the nu = 1 board:'
      % (kills, len(paths)))
print('the descent fires with r = 1 -- no opposed pair, no corridor,')
print('no template, no seal.')
