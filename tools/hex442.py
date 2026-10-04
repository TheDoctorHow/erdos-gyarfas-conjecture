"""hex442.py -- two decisive questions for the m = 3, t = 3 board.

 (1) ARC-FORCING PIGEONHOLE: in every hexagon (6 slots, spacing 4 on
     the 24-cycle), EVERY 2-colouring of the slots by leg parity has
     a monochromatic pair at arc 4 or arc 8.  (A colour class with
     all pairwise arcs = 12 sits inside one antipodal pair, size <=
     2; two such classes cover <= 4 < 6 slots.)  Machine: all 64
     colourings.  Consequence: the adversary can NEVER confine the
     even-harvest extras to arc 12 - the arc-12 SAT scenarios are
     unreachable.

 (2) DIAMETERS ALONE: is the 12-diameter board (no extras) already
     UNSAT under the full composite battery (subsets <= 5, all arc
     choices, mod 16), for the EVEN and for the ODD branch?  If yes,
     the extras are not even needed and the hexagon world dies
     without any placement casework.
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

# ---------------- (1) arc-forcing pigeonhole ----------------
HEXSLOTS = [0, 4, 8, 12, 16, 20]
def arc(p, q):
    d = abs(p - q) % L
    return min(d, L - d)
ok = 0
for bits in range(64):
    col = {HEXSLOTS[i]: (bits >> i) & 1 for i in range(6)}
    found = False
    for p, q in itertools.combinations(HEXSLOTS, 2):
        if col[p] == col[q] and arc(p, q) in (4, 8):
            found = True; break
    ok += found
print('(1) arc-forcing pigeonhole: %d/64 colourings have a'
      ' monochromatic pair at arc 4 or 8' % ok, flush=True)

# ---------------- (2) diameters alone, full battery ----------------
def gen_composites(ears, maxk=5):
    n = len(ears)
    for k in range(1, maxk + 1):
        for sub in itertools.combinations(range(n), k):
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

ears = [('d', u, (u + 12) % L) for u in range(12)]
t0 = time.time()
cons = sorted(set(gen_composites(ears, 5)))
print('(2) 12-diameter board, full battery maxk 5: %d constraints'
      ' (%.0fs to generate)' % (len(cons), time.time() - t0), flush=True)
# diameter shallowness: arc 12 -> T =/= 4 (mod 8)
def dom_mask(parity):
    m = 0
    for v in range(parity, 16, 2):
        if v % 8 != 4: m |= (1 << v)
    return m
for branch, parity in (('EVEN', 0), ('ODD', 1)):
    doms = [dom_mask(parity)] * 12
    t0 = time.time()
    res, sol, nd = solve_fc(doms, cons, 5 * 10 ** 7)
    extra = ('  T = %s' % sol) if res == 'SAT' else ''
    print('    %s branch: %-9s [%d nodes, %.1fs]%s'
          % (branch, res, nd, time.time() - t0, extra), flush=True)
    if res == 'SAT':
        okc = all((sum(sol[v] for v in vs) + kept) % 16 != 0
                  for vs, kept in cons)
        print('        witness verified: %s' % okc, flush=True)
