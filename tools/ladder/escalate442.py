"""escalate442.py -- the two-extra escalation for the ODD branch.

Context: the one-extra sweep (onextra442) kills every arc-4
placement in the odd-diameter branch, and the arc-forcing pigeonhole
(hex442) says every hexagon owns a same-parity pair at arc 4 or 8.
A hexagon offers ONLY arc-8 pairs exactly when its leg parities
alternate around the hexagon (the proper 2-colourings of a 6-cycle);
then its parity classes are {h, h+8, h+16} and {h+4, h+12, h+20},
each with three arc-8 pairs: six possible extra placements.

A surviving full board therefore has ALL FOUR hexagons alternating,
extras all at arc 8.  Every such board contains extras in hexagon 0
AND hexagon 1; its restriction to {12 diameters} + {those two
extras} is a sub-board, and sub-board UNSAT kills the full board.
So: decide all 6 x 6 slot combinations for hexagons (0, 1), odd
diameters, even extras, arc 8, full battery maxk 5.  If all 36 are
UNSAT, the odd branch of the hexagon world is EMPTY.
"""
import itertools, time

L = 24

def pc(x): return bin(x).count('1')

def arc_mask(a, b, choice):
    m = 0
    if choice == 0:
        j = a
        while j != b:
            m |= (1 << j); j = (j + 1) % L
    else:
        j = b
        while j != a:
            m |= (1 << j); j = (j + 1) % L
    return m

def gen_composites(ears, excl, maxk=5):
    n = len(ears)
    for k in range(1, maxk + 1):
        for sub in itertools.combinations(range(n), k):
            ok = True
            for idx in sub:
                if idx in excl and excl[idx] & set(sub): ok = False
            if not ok: continue
            feet = [(ears[i][1], ears[i][2]) for i in sub]
            for choices in itertools.product((0, 1), repeat=k):
                masks = [arc_mask(a, b, ch)
                         for (a, b), ch in zip(feet, choices)]
                xor = 0
                for m in masks: xor ^= m
                adj = {}
                def add(x, y):
                    adj.setdefault(x, []).append(y)
                    adj.setdefault(y, []).append(x)
                for j in range(L):
                    if not (xor >> j) & 1:
                        add(('c', j), ('c', (j + 1) % L))
                for ii, (a, b) in enumerate(feet):
                    add(('c', a), ('m', ii)); add(('m', ii), ('c', b))
                if any(len(v) != 2 for v in adj.values()): continue
                verts = list(adj); seen = {verts[0]}; st = [verts[0]]
                while st:
                    x = st.pop()
                    for y in adj[x]:
                        if y not in seen: seen.add(y); st.append(y)
                if len(seen) != len(verts): continue
                yield (sub, L - pc(xor))

def solve_fc(doms0, cons, budget):
    n = len(doms0)
    dom = list(doms0)
    val = [None] * n
    big = []
    for vs, kept in cons:
        if len(vs) == 1:
            v = vs[0]
            dom[v] &= ~(1 << ((-kept) % 16))
            if dom[v] == 0: return 'UNSAT', None, 0
        else:
            big.append((vs, kept))
    pairtab = [dict() for _ in range(n)]
    k3 = []
    for vs, kept in big:
        if len(vs) == 2:
            a, b = vs
            for (v, u) in ((a, b), (b, a)):
                tab = pairtab[v].setdefault(u, [0] * 16)
                for x in range(16):
                    tab[x] |= (1 << ((-kept - x) % 16))
        else:
            k3.append((vs, kept))
    by3 = [[] for _ in range(n)]
    for ci, (vs, kept) in enumerate(k3):
        for v in vs: by3[v].append(ci)
    cnt = [len(vs) for vs, _ in k3]
    ssum = [0] * len(k3)
    nodes = [0]
    def dfs():
        nodes[0] += 1
        if nodes[0] > budget: raise TimeoutError
        v = -1; bs = 99
        for u in range(n):
            if val[u] is None:
                s = pc(dom[u])
                if s < bs: bs = s; v = u
        if v < 0: return True
        cand = dom[v]; b = 0
        while cand:
            if cand & 1:
                val[v] = b
                trail = []; ok = True
                for u, tab in pairtab[v].items():
                    m = tab[b]
                    if val[u] is not None:
                        if (m >> val[u]) & 1: ok = False; break
                    else:
                        rem = dom[u] & m
                        if rem:
                            dom[u] &= ~rem; trail.append((u, rem))
                            if dom[u] == 0: ok = False; break
                touched = []
                if ok:
                    for ci in by3[v]:
                        cnt[ci] -= 1; ssum[ci] += b; touched.append(ci)
                        if cnt[ci] == 0:
                            if (ssum[ci] + k3[ci][1]) % 16 == 0:
                                ok = False; break
                        elif cnt[ci] == 1:
                            vs = k3[ci][0]
                            u = next(w for w in vs if val[w] is None)
                            forb = (-k3[ci][1] - ssum[ci]) % 16
                            if (dom[u] >> forb) & 1:
                                dom[u] &= ~(1 << forb)
                                trail.append((u, 1 << forb))
                                if dom[u] == 0: ok = False; break
                if ok and dfs(): return True
                for ci in touched:
                    cnt[ci] += 1; ssum[ci] -= b
                for (u, rem) in trail: dom[u] |= rem
                val[v] = None
            cand >>= 1; b += 1
        return False
    try:
        sat = dfs()
        return ('SAT' if sat else 'UNSAT'), (list(val) if sat else None), \
            nodes[0]
    except TimeoutError:
        return 'UNDECIDED', None, nodes[0]

def shallow_mask(arcd, parity):
    arcs = {arcd % 8, (L - arcd) % 8}
    m = 0
    for v in range(parity, 16, 2):
        if v % 8 not in arcs: m |= (1 << v)
    return m

t00 = time.time()
n_unsat = 0
alive = []
for s0 in (0, 4, 8, 12, 16, 20):
    row = []
    for s1 in (1, 5, 9, 13, 17, 21):
        ears = [('d', u, (u + 12) % L) for u in range(12)]
        ears.append(('e', s0, (s0 + 8) % L))
        ears.append(('e', s1, (s1 + 8) % L))
        excl = {12: {s0 % 12, (s0 + 8) % 12},
                13: {s1 % 12, (s1 + 8) % 12}}
        cons = sorted(set(gen_composites(ears, excl, 5)))
        doms = [shallow_mask(12, 1)] * 12 + \
               [shallow_mask(8, 0), shallow_mask(8, 0)]
        res, sol, nd = solve_fc(doms, cons, 5 * 10 ** 7)
        row.append('%s(%d)' % ('U' if res == 'UNSAT' else res[0], nd))
        if res == 'UNSAT':
            n_unsat += 1
        else:
            alive.append((s0, s1, res, sol))
    print('s0=%2d: %s   [%.0fs]' % (s0, '  '.join(row),
                                    time.time() - t00), flush=True)

print(flush=True)
print('two-extra (hexagons 0,1; arc 8; odd diameters, even extras):'
      ' %d/36 UNSAT' % n_unsat, flush=True)
if not alive:
    print('ALL 36 UNSAT: every alternating-parity full board dies via')
    print('its (hexagon 0, hexagon 1) sub-board.  Combined with the')
    print('one-extra arc-4 kills, the ODD branch of the m = 3, t = 3')
    print('hexagon world is EMPTY.')
else:
    print('%d surviving slot combinations (honest record):' % len(alive))
    for (s0, s1, res, sol) in alive:
        print('  s0=%2d s1=%2d: %s  T = %s' % (s0, s1, res, sol))
