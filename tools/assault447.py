"""assault447.py -- THE ASSAULT: the m = 3 endgame board at general t,
the constancy dichotomy, and the t-free forcing lemmas.

 (A) T-FREE FORCING.  At every t: the hexagon classes mod 2^{t-1}
     have SIX slots each (438); the slots split into two triangles
     with all internal arcs 2^t; every 2-colouring has a mono pair
     at arc 2^t exactly (the 442.3 pigeonhole, verbatim at slot
     level).  Machine: 64/64 colourings at t = 3, 4, 5.
 (B) CONSTANCY DICHOTOMY.
     (B1) ZIGZAG IMMUNITY: under the constant adversary (inner
          shifts == c0, outer == -c0), EVERY odd zigzag value is
          identically c0 - position-free.  Machine: all odd
          selections k <= 9 from a synthetic 12-braid.
     (B2) PAIR RIGIDITY: under constancy the pair value is
          2 c0 + 2(y_i - u_j); survival forces the cross-difference
          set {y_i - u_j : i < j} to avoid ONE class mod 2^t - the
          position rigidity that seeds the cascade.  Machine:
          formula verified on the synthetic braid.
 (C) THE (m = 3, t) ONE-EXTRA BOARDS, t = 3, 4, 5.
     L = 3 * 2^t; diameters (u, u + L/2), u = 0..L/2-1; one forced
     extra at (0, 2^t) (the triangle pigeonhole pair); exclusion of
     leg-sharing diameters; composites (single-cycle test) up to
     maxk; constraints  sum T + kept =/= 0  (mod 2^{t+1});
     shallowness  T =/= arc, L - arc  (mod 2^t).  Branches EVEN /
     ODD diameters, even extra.  t = 3 is the 442 revalidation;
     t = 4 and t = 5 are the uniformity experiment.
"""
import itertools, time, sys
sys.setrecursionlimit(100000)

# ---------------- (A) t-free triangle pigeonhole ----------------
print('=== (A) T-FREE TRIANGLE PIGEONHOLE ===', flush=True)
for t in (3, 4, 5):
    L = 3 * (1 << t); sp = 1 << (t - 1)
    slots = [k * sp for k in range(6)]
    def arc(p, q):
        d = abs(p - q) % L
        return min(d, L - d)
    tri1 = [slots[0], slots[2], slots[4]]
    tri2 = [slots[1], slots[3], slots[5]]
    arcs1 = sorted({arc(a, b) for a, b in itertools.combinations(tri1, 2)})
    ok = 0
    for bits in range(64):
        col = {slots[i]: (bits >> i) & 1 for i in range(6)}
        if any(col[p] == col[q] and arc(p, q) == (1 << t)
               for p, q in itertools.combinations(slots, 2)): ok += 1
    print('t=%d: L=%3d, slot spacing %2d, triangle arcs %s,'
          ' mono pair at arc 2^t: %d/64'
          % (t, L, sp, arcs1, ok), flush=True)

# ---------------- (B) constancy dichotomy ----------------
print(); print('=== (B) CONSTANCY DICHOTOMY ===', flush=True)
t = 4; L = 3 * (1 << t); c0 = 6
N = 12
us = [2 + 3 * i for i in range(N)]
ys = [L // 2 + 5 + 3 * i for i in range(N)]
# inner shift delta_i = c0 for all i; outer = c0 + 2 a_i - L == -c0
# requires 2 a_i == L - 2 c0; we only need the VALUE LAW, so define
# shifts symbolically:
inner = [c0] * N
outer = [-c0] * N
vals = set()
for k in range(1, 10, 2):
    for sel in itertools.combinations(range(N), k):
        v = sum(inner[i] if pos % 2 == 0 else outer[i]
                for pos, i in enumerate(sel))
        vals.add(v)
print('(B1) zigzag values over all odd selections k <= 9 of a'
      ' 12-braid under constancy: %s  (position-free: c0 = %d)'
      % (sorted(vals), c0), flush=True)
pair_ok = 0; tot = 0
for i in range(N):
    for j in range(i + 1, N):
        f_i = inner[i] + 2 * ys[i]          # f = delta + 2y
        g_j = inner[j] - 2 * us[j]          # g = delta - 2u
        v = f_i + g_j
        tot += 1
        if v == 2 * c0 + 2 * (ys[i] - us[j]): pair_ok += 1
print('(B2) pair value = 2 c0 + 2(y_i - u_j): %d/%d exact'
      ' => survival forbids ONE class of y_i - u_j mod 2^t'
      % (pair_ok, tot), flush=True)

# ---------------- (C) the general-t one-extra boards ----------------
print(); print('=== (C) THE (m=3, t) ONE-EXTRA BOARDS ===', flush=True)

def pcnt(x): return bin(x).count('1')

def board(t):
    L = 3 * (1 << t)
    ears = [('d', u, u + L // 2) for u in range(L // 2)]
    ears.append(('e', 0, 1 << t))
    excl = {len(ears) - 1: {0, (1 << t) % (L // 2)}}
    return L, ears, excl

def gen_cons(L, ears, excl, maxk):
    n = len(ears)
    seen = set()
    for k in range(1, maxk + 1):
        for sub in itertools.combinations(range(n), k):
            ok = True
            for idx in sub:
                if idx in excl and excl[idx] & set(sub): ok = False
            if not ok: continue
            feet = [(ears[i][1], ears[i][2]) for i in sub]
            for choices in itertools.product((0, 1), repeat=k):
                ced = 0
                for (a, b), ch in zip(feet, choices):
                    if ch == 0:
                        j = a
                        while j != b:
                            ced ^= (1 << j); j = (j + 1) % L
                    else:
                        j = b
                        while j != a:
                            ced ^= (1 << j); j = (j + 1) % L
                adj = {}
                def add(x, y):
                    adj.setdefault(x, []).append(y)
                    adj.setdefault(y, []).append(x)
                for j in range(L):
                    if not (ced >> j) & 1:
                        add(('c', j), ('c', (j + 1) % L))
                for ii, (a, b) in enumerate(feet):
                    add(('c', a), ('m', ii)); add(('m', ii), ('c', b))
                if any(len(v) != 2 for v in adj.values()): continue
                verts = list(adj); seenv = {verts[0]}; st = [verts[0]]
                while st:
                    x = st.pop()
                    for y in adj[x]:
                        if y not in seenv: seenv.add(y); st.append(y)
                if len(seenv) != len(verts): continue
                seen.add((sub, L - pcnt(ced)))
    return sorted(seen)

def dom_mask(L, t, a, b, parity):
    M = 1 << (t + 1)
    arcd = (b - a) % L
    forb = {arcd % (1 << t), (L - arcd) % (1 << t)}
    mask = 0
    for v in range(parity, M, 2):
        if v % (1 << t) not in forb: mask |= (1 << v)
    return mask

def solve_fc(doms0, cons, M, budget):
    n = len(doms0)
    dom = list(doms0)
    val = [None] * n
    big = []
    for vs, kept in cons:
        if len(vs) == 1:
            v = vs[0]
            dom[v] &= ~(1 << ((-kept) % M))
            if dom[v] == 0: return 'UNSAT', None, 0
        else:
            big.append((vs, kept))
    pairtab = [dict() for _ in range(n)]
    k3 = []
    for vs, kept in big:
        if len(vs) == 2:
            a, b = vs
            for (v, u) in ((a, b), (b, a)):
                tab = pairtab[v].setdefault(u, [0] * M)
                for x in range(M):
                    tab[x] |= (1 << ((-kept - x) % M))
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
        v = -1; bs = 10 ** 9
        for u in range(n):
            if val[u] is None:
                s = pcnt(dom[u])
                if s < bs: bs = s; v = u
        if v < 0: return True
        cand = dom[v]; b = 0
        while cand:
            if cand & 1:
                val[v] = b
                trail = []; ok = True
                for u, tab in pairtab[v].items():
                    mm = tab[b]
                    if val[u] is not None:
                        if (mm >> val[u]) & 1: ok = False; break
                    else:
                        rem = dom[u] & mm
                        if rem:
                            dom[u] &= ~rem; trail.append((u, rem))
                            if dom[u] == 0: ok = False; break
                touched = []
                if ok:
                    for ci in by3[v]:
                        cnt[ci] -= 1; ssum[ci] += b; touched.append(ci)
                        if cnt[ci] == 0:
                            if (ssum[ci] + k3[ci][1]) % M == 0:
                                ok = False; break
                        elif cnt[ci] == 1:
                            vs = k3[ci][0]
                            u = next(w for w in vs if val[w] is None)
                            forb = (-k3[ci][1] - ssum[ci]) % M
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

for t in (3, 4, 5):
    L, ears, excl = board(t)
    M = 1 << (t + 1)
    maxk = 3 if t >= 5 else 5
    t0 = time.time()
    cons = gen_cons(L, ears, excl, maxk)
    print('t=%d: L=%d, %d ears (one extra at (0,%d)), maxk=%d:'
          ' %d constraints (%.0fs)'
          % (t, L, len(ears), 1 << t, maxk, len(cons),
             time.time() - t0), flush=True)
    for branch, dpar in (('EVEN', 0), ('ODD', 1)):
        doms = [dom_mask(L, t, a, b, dpar) for (k, a, b) in ears[:-1]]
        doms.append(dom_mask(L, t, ears[-1][1], ears[-1][2], 0))
        t0 = time.time()
        res, sol, nd = solve_fc(doms, cons, M, 10 ** 8)
        extra = ('  T = %s' % sol) if res == 'SAT' else ''
        print('   %s branch: %-9s [%d nodes, %.1fs]%s'
              % (branch, res, nd, time.time() - t0, extra), flush=True)
        if res == 'SAT':
            okc = all((sum(sol[v] for v in vs) + kept) % M != 0
                      for vs, kept in cons)
            okd = all((doms[i] >> sol[i]) & 1 for i in range(len(sol)))
            print('       witness verified: constraints %s, domains %s'
                  % (okc, okd), flush=True)
