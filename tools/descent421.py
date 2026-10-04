"""descent421.py -- route (a): the Radical-Cut Descent.

Engine under test: in the knife-edge world (rank lambda = 2, plane <O_1,O_2>),
every monochromatic structure-cycle between consecutive same-side attachments
is RADICAL (meets O_1, O_2 evenly), hence a BICYCLE, hence an element of the
CUT SPACE. But a cycle that is a cut alternates sides along itself, which
confines every structure hanging at an odd-offset interior position of its
arc to attach ONLY inside the arc; the >= 3-attachment boundary law plus
coverage then force an infinite descent on strictly shrinking even intervals,
which bottoms out: SOME radical candidate is NOT a cut.

Machine corroboration (the theorem's conclusion, checked directly): in every
synthetic chordless monochrome board, at least one - in practice every -
same-side consecutive-attachment structure-cycle FAILS cut-space membership.
Cut space = orthogonal complement of the cycle space (edge dot product).

(A) one full handbuilt board (A = 15, O_2 = 7, two carriers, 5 windows);
(B) random chordless monochrome boards (D_3 + pockets, no O_2 needed for
    the engine), A in {11..21}: every candidate tested, zero boards where
    all candidates pass.
Boards are census-relaxed (mod-3 books not tuned): the engine under test
uses only monochromatism, coverage, >= 3 atts, and window laws.
"""
import itertools, random

def analyse(VER, EDG, cycles, tag, verbose=False):
    V = {v: i for i, v in enumerate(VER)}
    E = {}
    for e in EDG:
        E[frozenset(e)] = len(E)
    adj = {v: [] for v in VER}
    for u, w in EDG:
        adj[u].append(w); adj[w].append(u)
    bit = lambda u, w: 1 << E[frozenset((u, w))]
    # spanning tree -> fundamental cycles
    tree = set(); seen = {VER[0]}; stack = [VER[0]]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w); tree.add(frozenset((v, w))); stack.append(w)
    par = {VER[0]: None}; stack = [VER[0]]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if frozenset((v, w)) in tree and w not in par:
                par[w] = v; stack.append(w)
    def treepath(a, b):
        m = 0
        anc = set()
        x = a
        while x is not None: anc.add(x); x = par[x]
        y = b
        while y not in anc:
            m |= bit(y, par[y]); y = par[y]
        x = a
        while x != y:
            m |= bit(x, par[x]); x = par[x]
        return m
    fund = [(1 << E[e]) | treepath(*tuple(e)) for e in E if e not in tree]
    popc = lambda m: bin(m).count('1')
    incut = lambda z: all(popc(z & c) % 2 == 0 for c in fund)
    res = []
    for name, z in cycles:
        assert popc(z) % 2 == 0, name + ' not even'
        res.append((name, incut(z)))
    fails = [n for n, ok in res if not ok]
    if verbose:
        for name, ok in res:
            print('    %-28s %s' % (name, 'IN cut space' if ok
                                    else 'NOT a cut  -> KILL'))
    return len(res), len(fails)

def build_board(A, s1, WIN, pockets, carriers=None, B2=None):
    """O_1 cycle through b_1 with f = (s1, s1+1); D_3 caterpillar on WIN with
    b_3 spliced; pockets = list of att-lists; optional O_2 + carriers."""
    o = ['o%d' % i for i in range(A - 1)]
    VER = ['b1'] + o
    EDG = [('b1', 'o0'), ('o%d' % (A - 2), 'b1')] + \
          [(o[i], o[i + 1]) for i in range(A - 2)]
    # D_3 caterpillar: d_w - (b3 or z spacers) - d_w ...
    Ws = sorted(WIN)
    VER += ['d%d' % w for w in Ws] + ['b3']
    EDG += [('o%d' % w, 'd%d' % w) for w in Ws]
    spine = []
    for i in range(len(Ws) - 1):
        z = 'b3' if i == 0 else 'z%d' % i
        if i > 0: VER.append(z)
        spine.append(z)
        EDG += [('d%d' % Ws[i], z), (z, 'd%d' % Ws[i + 1])]
    cyc = []
    chi = lambda p: (p + (1 if p > s1 else 0)) % 2
    side = lambda p: 0 if p <= s1 else 1
    # D_3 consecutive same-side window cycles
    for i in range(len(Ws) - 1):
        w1, w2 = Ws[i], Ws[i + 1]
        if side(w1) != side(w2): continue
        z = {('o%d' % p, 'o%d' % (p + 1)) for p in range(w1, w2)}
        z |= {('o%d' % w1, 'd%d' % w1), ('o%d' % w2, 'd%d' % w2),
              ('d%d' % Ws[i], spine[i]), (spine[i], 'd%d' % Ws[i + 1])}
        cyc.append(('D3 windows (%d,%d)' % (w1, w2), z))
    # pockets: u-vertices take <= 2 atts + spine
    for k, atts in enumerate(pockets):
        atts = sorted(atts)
        us = []
        for j in range(0, len(atts), 2):
            u = 'P%du%d' % (k, j)
            VER.append(u); us.append((u, atts[j:j + 2]))
            EDG += [('o%d' % a, u) for a in atts[j:j + 2]]
        for j in range(len(us) - 1):
            z = 'P%dz%d' % (k, j)
            VER.append(z)
            EDG += [(us[j][0], z), (z, us[j + 1][0])]
        owner = {}
        for u, aa in us:
            for a in aa: owner[a] = u
        for i in range(len(atts) - 1):
            a1, a2 = atts[i], atts[i + 1]
            if side(a1) != side(a2) or chi(a1) != chi(a2): continue
            z = {('o%d' % p, 'o%d' % (p + 1)) for p in range(a1, a2)}
            z |= {('o%d' % a1, owner[a1]), ('o%d' % a2, owner[a2])}
            if owner[a1] != owner[a2]:
                i1 = us.index(next(x for x in us if x[0] == owner[a1]))
                i2 = us.index(next(x for x in us if x[0] == owner[a2]))
                for j in range(min(i1, i2), max(i1, i2)):
                    z |= {(us[j][0], 'P%dz%d' % (k, j)),
                          ('P%dz%d' % (k, j), us[j + 1][0])}
            cyc.append(('pocket%d atts (%d,%d)' % (k, a1, a2), z))
    if B2:
        p = ['p%d' % i for i in range(B2 - 2)]
        VER += ['b2'] + p
        EDG += [('b2', 'p0'), (p[-1], 'b2')] + \
               [(p[i], p[i + 1]) for i in range(B2 - 3)]
        for j, (foot, dest) in enumerate(carriers):
            c = 'c%d' % j
            VER.append(c)
            EDG += [('o%d' % foot, c), (c, 'p%d' % dest)]
    z2 = lambda s: {frozenset(e) for e in s}
    EDGf = []
    seenE = set()
    for e in EDG:
        if frozenset(e) not in seenE:
            seenE.add(frozenset(e)); EDGf.append(e)
    # convert cycle edge-sets to bitmasks later inside analyse via names
    return VER, EDGf, cyc

def to_masks(VER, EDG, cyc):
    E = {frozenset(e): i for i, e in enumerate(EDG)}
    out = []
    for name, z in cyc:
        m = 0
        for e in z: m |= 1 << E[frozenset(e)]
        out.append((name, m))
    return out

# ---------- (A) handbuilt full board: A = 15, s1 = 4, O_2 = 7, 2 carriers
VER, EDG, cyc = build_board(
    15, 4, WIN=[1, 3, 8, 10, 12],
    pockets=[[5, 7, 9], [0, 4, 11, 13]],
    carriers=[(2, 1), (6, 4)], B2=7)
print('(A) handbuilt board A=15 (windows {1,3,8,10,12}, pockets {5,7,9},'
      '{0,4,11,13}, carriers at 2,6):')
tot, fails = analyse(VER, EDG, to_masks(VER, EDG, cyc), 'A', verbose=True)
print('    %d/%d radical candidates are NOT cuts' % (fails, tot))
assert fails >= 1

# ---------- (B) random chordless monochrome boards
random.seed(421)
boards = all_fail = some_fail = 0
for trial in range(400):
    A = random.choice([11, 13, 15, 17, 19, 21])
    s1 = random.choice(range(2, A - 4, 2))
    chi = lambda p: (p + (1 if p > s1 else 0)) % 2
    pos = [p for p in range(A - 1)]
    c1 = [p for p in pos if chi(p) == 1]
    left1 = [p for p in c1 if p < s1]; right1 = [p for p in c1 if p > s1 + 1]
    if len(left1) + len(right1) < 3: continue
    if max(len(left1), len(right1)) < 2: continue
    kL = min(len(left1), random.choice([1, 2]))
    kR = 3 - kL if 3 - kL <= len(right1) else len(right1)
    if kL + kR < 3 or max(kL, kR) < 2:
        kL, kR = min(2, len(left1)), min(2, len(right1))
        if kL + kR < 3 or max(kL, kR) < 2: continue
    WIN = sorted(random.sample(left1, kL) + random.sample(right1, kR))
    rest = [p for p in pos if p not in WIN]
    byc = {0: [p for p in rest if chi(p) == 0],
           1: [p for p in rest if chi(p) == 1]}
    pockets = []
    bad = False
    for c in (0, 1):
        pool = byc[c][:]
        random.shuffle(pool)
        while pool:
            if len(pool) < 3: bad = True; break
            t = len(pool) if len(pool) < 6 else random.choice([3, 4, 5])
            pockets.append(sorted(pool[:t])); pool = pool[t:]
        if bad: break
    if bad: continue
    VER, EDG, cyc = build_board(A, s1, WIN, pockets)
    if not cyc: continue
    tot, fails = analyse(VER, EDG, to_masks(VER, EDG, cyc), 'B')
    boards += 1
    if fails == tot: all_fail += 1
    if fails >= 1: some_fail += 1
print()
print('(B) random boards: %d built; >= 1 candidate not a cut: %d/%d;'
      ' ALL candidates fail: %d/%d' % (boards, some_fail, boards,
                                       all_fail, boards))
assert some_fail == boards, 'DESCENT REFUTED on some board!'
print()
print('VERDICT: in every monochrome board tested, the knife-edge demand')
print('(all radical structure-cycles simultaneously cuts) is unsatisfiable.')
