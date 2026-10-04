"""diag442.py -- two-sided validation and diagnostics for the 442
decision.

 (a) SAT-side validation: EVEN arc 12 maxk 3 must come out SAT
     (matching 440's independent engine), witness fully verified.
 (b) Battery diagnostics: is the ODD branch already UNSAT at maxk 3
     (pairs + triples), or were quads/quintuples needed?
 (c) Independent consistency probe: 100000 random assignments from
     the ODD arc-4 domains -- every one must violate at least one
     constraint (consistent with UNSAT; sampling, not proof).
 (d) Completeness: ODD arc 12 maxk 5 (the scenario 440.2 shows moot;
     recorded anyway).
"""
import random, time
import importlib.util, sys, types

spec = importlib.util.spec_from_file_location(
    'cover442', 'tools/cover442.py')
# do not execute cover442's main body -- reimplement the small pieces
L = 24
import itertools

def build_board(extra_arc):
    ears = []
    for u in range(12):
        ears.append(('d', u, (u + 12) % L))
    for c in range(4):
        ears.append(('e', c, (c + extra_arc) % L))
    excl = {}
    for idx, (k, a, b) in enumerate(ears):
        if k == 'e':
            bad = {a % 12, b % 12}
            excl[idx] = {j for j in range(12) if j in bad}
    return ears, excl

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

def pc(x): return bin(x).count('1')

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
                yield (sub, choices, L - pc(xor), masks)

def shallowness_domain(kind, a, b, parity):
    arcd = (b - a) % L; arcs = {arcd % 8, (L - arcd) % 8}
    mask = 0
    for v in range(parity, 16, 2):
        if v % 8 not in arcs:
            mask |= (1 << v)
    return mask

def constraints_for(ears, excl, maxk):
    seen = set()
    for (sub, choices, kept, masks) in gen_composites(ears, excl, maxk):
        seen.add((sub, kept))
    return sorted(seen)

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

BUDGET = 5 * 10 ** 7
runs = [('(a) EVEN sat-side', 12, 3, {'d': 0, 'e': 0}),
        ('(b) ODD  battery',   4, 3, {'d': 1, 'e': 0}),
        ('(b) ODD  battery',   8, 3, {'d': 1, 'e': 0}),
        ('(d) ODD  arc 12',   12, 5, {'d': 1, 'e': 0})]
for label, extra_arc, maxk, pmap in runs:
    ears, excl = build_board(extra_arc)
    cons = constraints_for(ears, excl, maxk)
    doms = [shallowness_domain(k, a, b, pmap[k]) for (k, a, b) in ears]
    t0 = time.time()
    res, sol, nd = solve_fc(doms, cons, BUDGET)
    extra = ('  T = %s' % sol) if res == 'SAT' else ''
    print('%-18s arc %2d maxk %d: %4d constraints -> %-9s'
          ' [%d nodes, %.1fs]%s'
          % (label, extra_arc, maxk, len(cons), res, nd,
             time.time() - t0, extra), flush=True)
    if res == 'SAT':
        ok = all((sum(sol[v] for v in vs) + kept) % 16 != 0
                 for vs, kept in cons)
        okd = all((doms[i] >> sol[i]) & 1 for i in range(len(sol)))
        print('                   witness verified: constraints %s,'
              ' domains %s' % (ok, okd), flush=True)

# (c) random-assignment probe of the ODD arc-4 full battery
print(flush=True)
random.seed(442)
ears, excl = build_board(4)
cons = constraints_for(ears, excl, 5)
doms = [shallowness_domain(k, a, b, p)
        for (k, a, b), p in zip(ears, [1] * 12 + [0] * 4)]
vals = [[v for v in range(16) if (d >> v) & 1] for d in doms]
alive = 0
for _ in range(100000):
    x = [random.choice(vv) for vv in vals]
    if all((sum(x[v] for v in vs) + kept) % 16 != 0 for vs, kept in cons):
        alive += 1
print('(c) ODD arc 4 full battery: random probe 100000 assignments:'
      ' %d satisfy all constraints' % alive, flush=True)
