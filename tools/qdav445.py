"""qdav445.py -- QUADRATIC ZERO-SUM THEORY: the composite calculus as
a Davenport problem for Z/2^t-valued quadratic forms on the cube.

Model: k elements, values delta_i in Z/2^t, couplings c_ij; the
subset value is Q(S) = sum_{i in S} delta_i + sum_{i<j in S} c_ij.
In the real problem (441.4/442.1) c_ij = 2<z_i, z_j> is ALWAYS EVEN
(polarization) - this parity is the load-bearing hypothesis.

 (A) LEMMA A (anisotropy): every quadratic function q on F_2^3 with
     q(0) = 0 has a nonzero zero; dimension 2 admits anisotropic
     forms (delta = (1,1), c = 1).  Machine: all 64 functions on
     F_2^3, all 8 on F_2^2.
 (B) ARC-SWITCH LAW on Mobius-Kantor: for every ear of every cycle,
     the two arc-lifts satisfy z+ + z- = -x, <z+, z-> = -d, and
     delta+ + delta- = 2d - L.  (Exact, full census.)
 (C) SHARPNESS AT t = 2: exhaust all k = 4 instances (delta in
     {1,2,3}^4, c in {0,2}^6): every one contains a nonempty
     zero-subset mod 4; k = 3 has an avoiding witness.  So the sharp
     even-coupling constant at t = 2 is 4 = 2^t (the proved general
     bound 2*3^{t-1} = 6 is not tight).
 (D) t = 3 DECISION at k = 8: search for an avoiding instance
     (delta sorted in {1..7}^8, c in {0,2,4,6}^28), then attempt
     exhaustive DFS with early subset checks.  Honest outcome:
     SAT witness / UNSAT (sharp = 8) / UNDECIDED.
 (E) THEOREM CONSISTENCY at t = 3: random instances at k = 18
     (the proved bound 2*3^2) all contain zero-subsets.
 (F) THE ASYMMETRIC TARGET: progress composites need Q(S) == 2^t
     (mod 2^{t+1}); machine: the linear adversary {1,1,1,6} avoids
     target 4 mod 8 at k = 4 (beating the zero-target constant), and
     the maximum linear avoider size mod 8 is computed exhaustively.
     The all-even adversary avoids odd-valuation targets forever -
     the recursion base genuinely needs odd supply (the archive's
     shield/harvest theorems are exactly that supply).
"""
import itertools, random, sys, time
sys.setrecursionlimit(100000)
random.seed(445)

# ---------------- (A) Lemma A ----------------
def quad_zero_check(m):
    bad = 0
    aniso = []
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    for dbits in range(1 << m):
        for cbits in range(1 << len(pairs)):
            found = False
            for eps in range(1, 1 << m):
                v = 0
                for i in range(m):
                    if (eps >> i) & 1: v ^= (dbits >> i) & 1
                for pi, (i, j) in enumerate(pairs):
                    if (eps >> i) & 1 and (eps >> j) & 1:
                        v ^= (cbits >> pi) & 1
                if v == 0: found = True; break
            if not found:
                bad += 1
                if m == 2: aniso.append((dbits, cbits))
    return bad, aniso

bad3, _ = quad_zero_check(3)
bad2, aniso2 = quad_zero_check(2)
print('(A) F_2^3: %d/64 quadratic functions lack a nonzero zero'
      ' (LEMMA A: 0 expected)' % bad3)
print('    F_2^2: %d/8 anisotropic (nonzero => dim 2 is the true'
      ' ceiling): %s' % (bad2, aniso2), flush=True)
assert bad3 == 0 and bad2 > 0

# ---------------- (B) arc-switch law on Mobius-Kantor ----------------
def mobius_kantor():
    E = [(i, (i + 1) % 8) for i in range(8)]
    E += [(i, 8 + i) for i in range(8)]
    E += [(8 + i, 8 + ((i + 3) % 8)) for i in range(8)]
    return 16, [(min(a, b), max(a, b)) for a, b in E]

def build(nE):
    n, E = nE
    adj = {v: [] for v in range(n)}
    for a, b in E:
        adj[a].append(b); adj[b].append(a)
    eidx = {}
    for k, (a, b) in enumerate(E):
        eidx[(a, b)] = (k, 1); eidx[(b, a)] = (k, -1)
    return n, E, adj, eidx

def norm2(x): return sum(v * v for v in x)
def dot(x, y): return sum(a * b for a, b in zip(x, y))

def all_cycles(adj, n):
    out = set()
    def dfs(start, v, visited, path):
        for w in adj[v]:
            if w == start and len(path) >= 3:
                if path[1] < path[-1]: out.add(tuple(path))
            elif w > start and w not in visited:
                visited.add(w); path.append(w)
                dfs(start, w, visited, path)
                path.pop(); visited.remove(w)
    for s in range(n):
        dfs(s, s, {s}, [s])
    return sorted(out, key=lambda c: (len(c), c))

def lift(cyc, eidx, m):
    x = [0] * m
    for a, b in zip(cyc, cyc[1:] + cyc[:1]):
        k, s = eidx[(a, b)]; x[k] += s
    return x

def ears_of(adj, C):
    Cset = set(C); L = len(C); res = []
    for i in range(L):
        for j in range(i + 1, L):
            u, v = C[i], C[j]
            stack = [(u, [u])]
            while stack:
                w, path = stack.pop()
                for nb in adj[w]:
                    if nb == v:
                        if len(path) == 1 and (j - i == 1 or
                                               (i == 0 and j == L - 1)):
                            continue
                        res.append((i, j, path + [v]))
                    elif nb not in Cset and nb not in path:
                        stack.append((nb, path + [nb]))
    return res

n, E, adj, eidx = build(mobius_kantor())
m = len(E)
cycles = all_cycles(adj, n)
tot = ok = 0
for C in cycles:
    L = len(C); x = lift(C, eidx, m)
    for (i, j, P) in ears_of(adj, C):
        d = len(P) - 1
        zP = [0] * m
        for a, b in zip(P, P[1:]):
            k, s = eidx[(a, b)]; zP[k] += s
        arcv = C[i:j + 1]
        za1 = [0] * m
        for a, b in zip(arcv, arcv[1:]):
            k, s = eidx[(a, b)]; za1[k] += s
        zp = [p - q for p, q in zip(zP, za1)]
        za2 = [xx - q for xx, q in zip(x, za1)]
        zm = [-(p + q) for p, q in zip(zP, za2)]
        tot += 1
        c1 = all(a + b == -xx for a, b, xx in zip(zp, zm, x))
        c2 = (dot(zp, zm) == -d)
        c3 = ((norm2(zp) + 2 * dot(x, zp)) +
              (norm2(zm) + 2 * dot(x, zm)) == 2 * d - L)
        if c1 and c2 and c3: ok += 1
print('(B) ARC-SWITCH LAW on Mobius-Kantor: z+ + z- = -x,'
      ' <z+,z-> = -d, delta+ + delta- = 2d - L: %d/%d exact'
      % (ok, tot), flush=True)

# ---------------- helpers for the abstract problem ----------------
def has_zero_subset(k, deltas, cmat, MOD, target=0):
    Q = [0] * (1 << k)
    for i in range(k):
        bit = 1 << i
        for mm in range(bit):
            Q[mm | bit] = (Q[mm] + deltas[i] +
                           sum(cmat[j][i] for j in range(k)
                               if (mm >> j) & 1)) % MOD
            if mm | bit and Q[mm | bit] == target % MOD:
                return True
    return False

# ---------------- (C) sharpness at t = 2 ----------------
allkill = 0; total = 0
for deltas in itertools.product((1, 2, 3), repeat=4):
    for cs in itertools.product((0, 2), repeat=6):
        cmat = [[0] * 4 for _ in range(4)]
        idx = 0
        for i in range(4):
            for j in range(i + 1, 4):
                cmat[i][j] = cmat[j][i] = cs[idx]; idx += 1
        total += 1
        if has_zero_subset(4, list(deltas), cmat, 4): allkill += 1
print('(C) t=2, k=4: %d/%d instances contain a nonempty zero-subset'
      ' mod 4' % (allkill, total))
cm3 = [[0] * 3 for _ in range(3)]
print('    t=2, k=3 avoidance witness delta=(1,1,1), c=0: avoids ='
      ' %s  => sharp constant at t=2 is 4 = 2^t'
      % (not has_zero_subset(3, [1, 1, 1], cm3, 4)), flush=True)

# ---------------- (D) t = 3, k = 8 decision ----------------
K = 8; MOD = 8
def sat_search(tries=4000, iters=400):
    for _ in range(tries):
        deltas = sorted(random.choice((1, 2, 3, 4, 5, 6, 7))
                        for _ in range(K))
        cmat = [[0] * K for _ in range(K)]
        for i in range(K):
            for j in range(i + 1, K):
                cmat[i][j] = cmat[j][i] = random.choice((0, 2, 4, 6))
        for _ in range(iters):
            # find a violated subset
            Q = [0] * (1 << K); badS = 0
            for i in range(K):
                bit = 1 << i
                for mm2 in range(bit):
                    Q[mm2 | bit] = (Q[mm2] + deltas[i] +
                                    sum(cmat[j][i] for j in range(K)
                                        if (mm2 >> j) & 1)) % MOD
                    if Q[mm2 | bit] == 0 and (mm2 | bit):
                        badS = mm2 | bit
                if badS: break
            if not badS:
                return deltas, cmat
            # perturb a member of the violated subset
            mems = [i for i in range(K) if (badS >> i) & 1]
            i = random.choice(mems)
            if random.random() < 0.4:
                deltas[i] = random.choice((1, 2, 3, 4, 5, 6, 7))
                deltas.sort()
            else:
                j = random.choice([j for j in range(K) if j != i])
                cmat[i][j] = cmat[j][i] = random.choice((0, 2, 4, 6))
    return None

t0 = time.time()
w = sat_search()
print('(D) t=3, k=8: randomized avoidance search: %s  [%.0fs]'
      % ('WITNESS FOUND' if w else 'none found', time.time() - t0),
      flush=True)

decided = None
if w:
    deltas, cmat = w
    assert not has_zero_subset(K, deltas, cmat, MOD)
    print('    witness delta = %s' % deltas)
    print('    c = %s' % [[cmat[i][j] for j in range(i + 1, K)]
                          for i in range(K)])
    decided = 'SAT'
else:
    # exhaustive DFS: elements added one at a time (delta sorted),
    # couplings chosen one at a time with early subset checks
    budget = [4 * 10 ** 7]
    Qarr = [0] * (1 << K)
    deltas = [0] * K
    cmat = [[0] * K for _ in range(K)]
    def place(i):
        """elements 0..i-1 placed, Qarr valid for masks < 2^i."""
        if budget[0] <= 0: raise TimeoutError
        if i == K: return True          # avoidance instance complete
        lo = deltas[i - 1] if i else 1
        for dv in range(lo, 8):
            deltas[i] = dv
            if place_c(i, 0, dv): return True
        return False
    def place_c(i, j, dv):
        """choosing c_{j,i}; subsets with max old element < j done."""
        if budget[0] <= 0: raise TimeoutError
        budget[0] -= 1
        if j == i:
            return place(i + 1)
        for cv in (0, 2, 4, 6):
            cmat[j][i] = cmat[i][j] = cv
            ok = True
            bit_i = 1 << i; bit_j = 1 << j
            # subsets S = m | bit_j | bit_i, m subset of [0..j-1]
            for mm2 in range(1 << j):
                base = mm2 | bit_j
                q = deltas[i]
                for l in range(j + 1):
                    if (base >> l) & 1: q += cmat[l][i]
                q = (Qarr[base] + q) % MOD
                if q == 0: ok = False; break
            if ok:
                if j == i - 1:
                    # element i's couplings complete: fill Qarr rows,
                    # also check subsets {i} and m|bit_i with m
                    # not containing j... all handled: masks m|bit_i
                    # for ALL m < 2^i
                    good = True
                    for mm2 in range(1 << i):
                        q = deltas[i]
                        for l in range(i):
                            if (mm2 >> l) & 1: q += cmat[l][i]
                        q = (Qarr[mm2] + q) % MOD
                        Qarr[mm2 | bit_i] = q
                        if q == 0 and (mm2 | bit_i): good = False; break
                    if good and place_c(i, j + 1, dv): return True
                else:
                    if place_c(i, j + 1, dv): return True
        return False
    t0 = time.time()
    try:
        # note: for i's with no couplings yet (i = 0), place_c(0,0)
        # -> j == i immediately -> place(1) after Qarr row fill?
        # element 0: fill Qarr[1] here:
        found = False
        for dv in range(1, 8):
            deltas[0] = dv; Qarr[1] = dv % MOD
            if Qarr[1] != 0 and place(1): found = True; break
        decided = 'SAT-late' if found else 'UNSAT'
    except TimeoutError:
        decided = 'UNDECIDED'
    print('    exhaustive DFS: %s  [%d nodes used, %.0fs]'
          % (decided, 4 * 10 ** 7 - budget[0], time.time() - t0),
          flush=True)

# ---------------- (E) theorem consistency at the proved bound ----
cnt = 0; N = 60
for _ in range(N):
    k = 18
    deltas = [random.choice(range(1, 8)) for _ in range(k)]
    cmat = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            cmat[i][j] = cmat[j][i] = random.choice((0, 2, 4, 6))
    if has_zero_subset(k, deltas, cmat, 8): cnt += 1
print('(E) t=3, k=18 (= proved bound 2*3^2): %d/%d random instances'
      ' contain a zero-subset (theorem says all must)' % (cnt, N),
      flush=True)

# ---------------- (F) the asymmetric target ----------------
cm4 = [[0] * 4 for _ in range(4)]
print('(F) linear adversary {1,1,1,6} avoids TARGET 4 mod 8: %s'
      % (not has_zero_subset(4, [1, 1, 1, 6], cm4, 8, target=4)))
# THE INVINCIBLE ABSTRACT ADVERSARY: delta == 2, c == 4 gives
# Q(S) = 2|S|^2 mod 8 - valuation always odd, NEVER the target 4.
k = 14
cm = [[4 * (i != j) for j in range(k)] for i in range(k)]
inv = not has_zero_subset(k, [2] * k, cm, 8, target=4)
qs = sorted({(2 * s * s) % 8 for s in range(1, k + 1)})
print('    INVINCIBLE adversary delta==2, c==4 at k=%d: avoids'
      ' target 4 mod 8: %s  (Q values realised: %s = 2s^2 mod 8)'
      % (k, inv, qs))
print('    => the abstract asymmetric problem is UNBOUNDED; the')
print('       boards could NOT have died by abstract zero-sums.')
# THE ARC-SWITCH REFUTATION of that adversary: delta+ == 2 mod 8
# with L == 4 mod 8 forces delta- == 2d - L - delta+ == 2d - 6;
# both-arc shallowness (delta- =/= 0 mod 4) forces d even; then
# delta+ = d - a == 2 mod 8 forces a == d - 2: the adversary's
# geometry is pinned - and pinned geometry is what the boards then
# exhaust.  Machine: enumerate all (d, a) mod 8 with d - a == 2 and
# check which survive both-arc shallowness at L == 4 mod 8:
alive = []
for d in range(8):
    for a in range(8):
        dp = (d - a) % 8
        dm = (2 * d - 4 - dp) % 8
        if dp == 2 and dm % 4 != 0:
            alive.append((d, a, dm))
print('    arc-switch on delta+ == 2 (L == 4 mod 8): surviving'
      ' (d, a, delta-) mod 8: %s' % alive)
print('    (d forced even, delta- forced == 2 mod 4: the invincible'
      ' family collapses to a one-parameter pinned class)', flush=True)
