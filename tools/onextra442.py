"""onextra442.py -- THE ONE-EXTRA REDUCTION for the m = 3, t = 3
hexagon world.

Logic: (i) the 12 diameters are forced (439); (ii) by the arc-forcing
pigeonhole (hex442: 64/64), every hexagon contains a same-leg-parity
slot pair at arc 4 or arc 8, and one banked even-harvest extra needs
only a tree path, so EVERY configuration's board contains the 12
diameters plus at least one extra at some (s, s+a), a in {4, 8};
(iii) the sub-board's composites are composites of any fuller board,
so sub-board UNSAT kills every superset.

Therefore: if the 13-ear board {12 diameters} + {extra at (s, s+a)}
is UNSAT for ALL s in 0..23 and a in {4, 8}, in a given parity
branch, that branch of the hexagon world is EMPTY - mixed scenarios,
all placements, no casework.

Extra parity: the banked extra joins two SAME-parity legs, so its
length T = l + l' + path has the parity of its connecting path; both
parities are run (extra even AND extra odd) so no parity assumption
on the extra survives unexamined.  48 placements x 2 diameter
branches x 2 extra parities = 192 exhaustive decisions.
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
worst = {}
survivors = []
for dpar in (0, 1):
    for epar in (0, 1):
        results = {}
        for a in (4, 8):
            for s in range(24):
                ears = [('d', u, (u + 12) % L) for u in range(12)]
                ears.append(('e', s, (s + a) % L))
                excl = {12: {s % 12, (s + a) % 12}}
                cons = sorted(set(gen_composites(ears, excl, 5)))
                doms = [shallow_mask(12, dpar)] * 12 + \
                       [shallow_mask(a, epar)]
                res, sol, nd = solve_fc(doms, cons, 5 * 10 ** 7)
                results[(a, s)] = (res, nd, len(cons))
                if res != 'UNSAT':
                    survivors.append((dpar, epar, a, s, res, sol))
        n_unsat = sum(1 for r in results.values() if r[0] == 'UNSAT')
        mx = max(r[1] for r in results.values())
        print('diam parity %d, extra parity %d: %d/48 placements UNSAT'
              '  (max nodes %d)  [%.0fs]'
              % (dpar, epar, n_unsat, mx, time.time() - t00), flush=True)

print(flush=True)
if not survivors:
    print('ALL 192 ONE-EXTRA BOARDS UNSAT: by the arc-forcing')
    print('pigeonhole and the sub-board principle, the m = 3, t = 3')
    print('hexagon world is EMPTY - both diameter parities, every')
    print('extra placement, every mixed scenario.')
else:
    print('%d SURVIVING one-extra boards (honest record):' % len(survivors))
    for (dpar, epar, a, s, res, sol) in survivors[:20]:
        print('  dpar %d epar %d arc %d slot %2d: %s  T = %s'
              % (dpar, epar, a, s, res, sol))
