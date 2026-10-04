"""load440.py -- the fully loaded diameter+extras board at t=3, m=3.

Board: C = 24-cycle; 12 diameter ears (u, u+12); 4 extras, one per
hexagon (same-class pairs at arc a_E in {4, 8, 12}), each extra
sharing legs with exactly two diameters (excluded from common
subsets). Composites: every subset of pairwise interior-disjoint
ears up to size 5, every arc-choice vector; single-cycle composites
(topological XOR test with dummy paths) give constraints
    sum of used T_i  =/=  -const  (mod 16)
(progress avoidance: composite length == 0 mod 16 would exceed
v_2(24) = 3). Per-ear shallowness: T =/= 4 (mod 8) for diameters
(arcs 12/12) and for extras with arcs {4,20},{8,16},{12,12}: all
reduce to T =/= arc, -arc (mod 8).

CSP: T_i mod 16 over EVEN or ODD domains; DFS with incremental
constraint checks. Output per branch and extra-arc scenario:
SAT witness or UNSAT (exhaustive).
"""
import itertools, sys
sys.setrecursionlimit(10000)

L = 24
def build_board(extra_arc):
    ears = []
    for u in range(12):
        ears.append(('d', u, (u + 12) % L))
    for c in range(4):
        ears.append(('e', c, (c + extra_arc) % L))
    # exclusion: extra c shares legs with diameters at feet c and
    # c+extra_arc (their leg-owning diameters: u = c mod 12 and
    # u = (c+extra_arc) mod 12)
    excl = {}
    for idx, (k, a, b) in enumerate(ears):
        if k == 'e':
            bad = {a % 12, b % 12}
            excl[idx] = {j for j in range(12) if j in bad}
    return ears, excl

def composites(ears, excl, maxk=5):
    """Yield (used_indices, const) for single-cycle composites."""
    n = len(ears)
    # edge universe: C edges 0..23 (i -> i+1), plus each ear a dummy
    # 2-edge path with a private middle vertex
    out = []
    for k in range(1, maxk + 1):
        for sub in itertools.combinations(range(n), k):
            ok = True
            for idx in sub:
                if idx in excl and excl[idx] & set(sub): ok = False
            if not ok: continue
            feet = [(ears[i][1], ears[i][2]) for i in sub]
            for choices in itertools.product((0, 1), repeat=k):
                # arc for ear i: choice 0 = arc from a forward to b,
                # 1 = the complement
                cedges = [0] * L
                for (a, b), ch in zip(feet, choices):
                    if ch == 0:
                        j = a
                        while j != b:
                            cedges[j] ^= 1; j = (j + 1) % L
                    else:
                        j = b
                        while j != a:
                            cedges[j] ^= 1; j = (j + 1) % L
                # composite edge set: C xor (arcs) plus ear paths
                # C edge j present iff 1 xor cedges[j]
                adj = {}
                def add(x, y):
                    adj.setdefault(x, []).append(y)
                    adj.setdefault(y, []).append(x)
                for j in range(L):
                    if 1 ^ cedges[j]:
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
                # length = sum T_i + (kept C edges) ... each ear path
                # contributes T_i (dummy 2 edges stand for T_i)
                kept = sum(1 ^ cedges[j] for j in range(L))
                out.append((sub, kept))
    # deduplicate constraint patterns
    return sorted(set(out))

def shallowness_domain(kind, a, b, parity):
    arcd = (b - a) % L; arcs = {arcd % 8, (L - arcd) % 8}
    dom = []
    for v in range(parity, 16, 2):
        if v % 8 not in arcs:
            dom.append(v)
    return dom

def solve(ears, excl, cons, parity_map, budget=3 * 10 ** 6):
    import random as rnd
    n = len(ears)
    doms = [shallowness_domain(k, a, b, parity_map[k])
            for (k, a, b) in ears]
    # (a) randomized search for a witness
    byvar = [[] for _ in range(n)]
    for sub, kept in cons:
        for j in sub: byvar[j].append((sub, kept))
    rnd.seed(440)
    for _ in range(150):
        sol = [rnd.choice(doms[i]) for i in range(n)]
        for _ in range(80):
            bad = [(sub, kept) for sub, kept in cons
                   if (sum(sol[j] for j in sub) + kept) % 16 == 0]
            if not bad:
                return 'SAT', sol
            sub, kept = rnd.choice(bad)
            j = rnd.choice(sub)
            sol[j] = rnd.choice(doms[j])
    # (b) budgeted exhaustive DFS
    bymax = {}
    for sub, kept in cons:
        bymax.setdefault(max(sub), []).append((sub, kept))
    sol = [None] * n
    nodes = [0]
    def dfs(i):
        if nodes[0] > budget: raise TimeoutError
        nodes[0] += 1
        if i == n: return True
        for v in doms[i]:
            sol[i] = v
            ok = True
            for sub, kept in bymax.get(i, []):
                tot = sum(sol[j] for j in sub) + kept
                if tot % 16 == 0: ok = False; break
            if ok and dfs(i + 1): return True
        sol[i] = None
        return False
    try:
        sat = dfs(0)
        return ('SAT' if sat else 'UNSAT'), sol
    except TimeoutError:
        return 'UNDECIDED', None

import time
for extra_arc in (4, 8, 12):
    ears, excl = build_board(extra_arc)
    for maxk in (3, 5):
        t0 = time.time()
        cons = composites(ears, excl, maxk=maxk)
        print('extra arc %2d, maxk=%d: %d constraints (%.0fs)'
              % (extra_arc, maxk, len(cons), time.time() - t0),
              flush=True)
        alive = False
        for branch, pmap in (('EVEN', {'d': 0, 'e': 0}),
                             ('ODD',  {'d': 1, 'e': 0})):
            res, sol = solve(ears, excl, cons, pmap)
            alive |= (res != 'UNSAT')
            print('  %s branch: %s%s'
                  % (branch, res,
                     ('  T = ' + str(sol)) if res == 'SAT' else ''),
                  flush=True)
        if not alive:
            print('  both branches UNSAT at maxk=%d' % maxk)
            break
