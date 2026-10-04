"""harvest438.py -- the 4/3 gap: the exact even-harvest function h(r)
and the interleaving dividend.

 (1) h(r) = worst case, over subcubic trees with r leaves and over
     2-colourings of the leaves, of the maximum number of pairwise
     vertex-disjoint same-colour leaf-pairs (tree paths). Each such
     pair is an even-shift ear (437.3). Exhaustive-ish: all sampled
     subcubic trees (caterpillars + random) for r = 3..8, all
     colourings, exact matching by brute force.
 (2) worst-case total harvest: minimise sum h(r_H) over compositions
     of L into parts >= 3 - the composition adversary's best split.
 (3) INTERLEAVING DIVIDEND: even extras as odd-scale deviants in the
     zigzag AP: with R alternation runs the reachable j-range is
     ~R/2; progress needs j* in range, j* up to 2^(t+1) - 1 for odd
     difference: required runs ~2^(t+2) worst case; available runs
     <= 2 * harvest. Table of the resulting factor at m = 3.
"""
import itertools, random

random.seed(438)

def rand_subcubic_tree(r):
    adj = {0: [1], 1: [0]}; nxt = 2
    while True:
        leaves = [v for v in adj if len(adj[v]) == 1]
        if len(leaves) == r: return adj
        if len(leaves) > r: return None
        cands = [v for v in adj if len(adj[v]) < 3]
        v = random.choice(cands)
        adj[v].append(nxt); adj[nxt] = [v]; nxt += 1

def caterpillar(r, leglens):
    adj = {}
    def add(a, b):
        adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    spine = ['s%d' % i for i in range(r)]
    for i in range(r - 1): add(spine[i], spine[i + 1])
    for i in range(r):
        prev = spine[i]
        for j in range(leglens[i % len(leglens)]):
            v = 'l%d_%d' % (i, j); add(prev, v); prev = v
    return adj

def tree_path(adj, a, b):
    par = {a: None}; st = [a]
    while st:
        v = st.pop()
        if v == b: break
        for w in adj[v]:
            if w not in par: par[w] = v; st.append(w)
    path = [b]
    while par[path[-1]] is not None: path.append(par[path[-1]])
    return set(path)

def max_disjoint_same(adj, colours):
    leaves = sorted([v for v in adj if len(adj[v]) == 1])
    pairs = [(a, b) for a, b in itertools.combinations(leaves, 2)
             if colours[a] == colours[b]]
    paths = {p: tree_path(adj, *p) for p in pairs}
    best = 0
    def rec(idx, used, count):
        nonlocal best
        best = max(best, count)
        if idx >= len(pairs): return
        if count + (len(pairs) - idx) <= best: return
        p = pairs[idx]
        if not (paths[p] & used):
            rec(idx + 1, used | paths[p], count + 1)
        rec(idx + 1, used, count)
    rec(0, set(), 0)
    return best

print('(1) the even-harvest function h(r):')
H = {}
for r in range(3, 8):
    worst = 99
    trees = [caterpillar(r, [1]), caterpillar(r, [1, 2]),
             caterpillar(r, [2])]
    seen = 0
    while seen < 40:
        T = rand_subcubic_tree(r)
        if T: trees.append(T); seen += 1
    for T in trees:
        leaves = sorted([v for v in T if len(T[v]) == 1])
        for bits in range(1 << (r - 1)):        # colourings up to flip
            colours = {leaves[i]: (bits >> i) & 1 if i < r - 1 else 0
                       for i in range(r)}
            worst = min(worst, max_disjoint_same(T, colours))
        if worst == 0: break
    H[r] = worst
    print('    h(%d) = %d' % (r, worst))

# (2) worst composition
print('(2) worst-case total harvest over compositions of L (parts 3..7):')
for t in (3, 4, 5):
    L = 3 * (1 << t)
    # min sum H over compositions: DP
    INF = 10 ** 9
    dp = [INF] * (L + 1); dp[0] = 0
    for s in range(1, L + 1):
        for r, hr in H.items():
            if s >= r: dp[s] = min(dp[s], dp[s - r] + hr)
    f1 = (1 << t) - 1
    print('    t=%d: L=%3d  worst harvest = %3d  vs f(1) = %3d '
          ' (ratio %.2f)' % (t, L, dp[L], f1, dp[L] / f1))

# (3) interleaving dividend
print('(3) interleaving dividend at m = 3: runs available vs required:')
for t in (3, 4, 5):
    L = 3 * (1 << t)
    INF = 10 ** 9
    dp = [INF] * (L + 1); dp[0] = 0
    for s in range(1, L + 1):
        for r, hr in H.items():
            if s >= r: dp[s] = min(dp[s], dp[s - r] + hr)
    runs_avail = 2 * dp[L]
    runs_need = 1 << (t + 2)
    print('    t=%d: runs available <= %3d, worst-case required %3d'
          ' (ratio %.2f)' % (t, runs_avail, runs_need,
                             runs_avail / runs_need))
print()
print('The gap is NOT 4/3: the worst-case harvest (r = 5 and r = 7'
      ' bridges) lands at L/5..2L/7, and the interleaving requirement')
print('is a further factor away. RECORDED HONESTLY: 437.4 was'
      ' optimistic; the razor is thicker than one harvest improvement.')
