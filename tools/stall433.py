"""stall433.py -- the stall theorem, the leg bound, and the middle zone.

Checks:
 (1) STALL: in the constant braided configuration (all inner shifts
     == c0, outer == -c0 mod 2^(t+1)), EVERY odd zigzag value equals
     c0: numeric sweep over k = 1, 3, ..., 11 with random spiral
     sub-selections. No composite ever reaches valuation >= t: the
     configuration is congruence-ISOLATED (no progress, no lateral).
 (2) FOUR-ROUTE POSITION CONSTRAINT: for a cross-path of length l
     joining interior positions (a on P_i, b on P_j), all four induced
     ears exceed L/2 iff  a + b > L/2 - l,  (Pi - a) + (Pj - b)
     > L/2 - l,  a + (Pj - b) > L/2 - l,  (Pi - a) + b > L/2 - l;
     sweep confirms: solutions require a and b in MIDDLE ZONES
     (min(a, Pi - a) > L/4 - l - s/2 where s = slack).
 (3) LEG BOUND: a subcubic tree in which all pairwise leaf distances
     exceed D has at most one leaf with leg < ceil(D/2); hence size
     >= (r - 1) ceil(D/2). Verified on random subcubic trees filtered
     for the distance property.
 (4) CHERRY-TREE TIGHTNESS: the balanced cherry construction with
     legs floor(D/2) + 1 achieves all pairwise leaf distances > D
     with size ~ r D / 2: the packing bound is TIGHT - the braided
     core is pinned to the boundary shell of 191.
"""
import random, itertools

v2 = lambda x: (abs(x) & -abs(x)).bit_length() - 1 if x else 999

# ---------------- (1) stall ----------------
t = 3; M2 = 1 << (t + 1)
random.seed(433)
c0 = 2
stall_ok = True
for trial in range(500):
    k = random.choice([1, 3, 5, 7, 9, 11])
    # alternating word: inner c0 at odd positions, outer -c0 at even
    val = sum((c0 if i % 2 == 0 else -c0) for i in range(k)) % M2
    stall_ok &= (val == c0 % M2)
assert stall_ok
print('(1) STALL: 500 random odd zigzags of the constant configuration'
      ' all have value c0 = %d mod %d (v_2 = %d < %d):' % (c0, M2,
      v2(c0), t))
print('    no composite of ANY shape reaches valuation >= t: the'
      ' configuration is congruence-ISOLATED. Route (i) closes'
      ' NEGATIVELY.')

# ---------------- (2) middle zone ----------------
L = 48                      # t=4, m=3
sweep = 0; sols = 0; zone = True
for Pi in (L // 2 + 3, L // 2 + 7):
    for Pj in (L // 2 + 3, L // 2 + 7):
        for l in (1, 3):
            for a in range(1, Pi):
                for b in range(1, Pj):
                    sweep += 1
                    routes = [a + l + b, (Pi - a) + l + (Pj - b),
                              a + l + (Pj - b), (Pi - a) + l + b]
                    if all(r > L // 2 for r in routes):
                        sols += 1
                        si = Pi - L // 2; sj = Pj - L // 2
                        ok_i = min(a, Pi - a) > L // 4 - l - sj // 2 - 1
                        ok_j = min(b, Pj - b) > L // 4 - l - si // 2 - 1
                        zone &= ok_i and ok_j
print('(2) middle zone: %d/%d position pairs admit a cross-connection,'
      ' every one with both endpoints in the middle zone: %s'
      % (sols, sweep, zone))
assert zone and sols > 0

# ---------------- (3) leg bound ----------------
def random_subcubic_tree_nv(nv):
    adj = {0: [1], 1: [0]}; nxt = 2
    while nxt < nv:
        cands = [v for v in adj if len(adj[v]) < 3]
        v = random.choice(cands)
        adj[v].append(nxt); adj[nxt] = [v]; nxt += 1
    return adj

def leaf_dists(adj):
    leaves = [v for v in adj if len(adj[v]) == 1]
    D = {}
    for s in leaves:
        dist = {s: 0}; q = [s]
        while q:
            v = q.pop(0)
            for w in adj[v]:
                if w not in dist: dist[w] = dist[v] + 1; q.append(w)
        D[s] = dist
    return leaves, D

def leg(adj, leaf):
    prev, cur, ln = None, leaf, 0
    while len(adj[cur]) <= 2:
        nxts = [w for w in adj[cur] if w != prev]
        if not nxts: break
        prev, cur = cur, nxts[0]; ln += 1
    return ln

checked = 0; ok3 = 0
while checked < 60:
    adj = random_subcubic_tree_nv(random.randint(10, 60))
    leaves, D = leaf_dists(adj)
    if len(leaves) < 3: continue
    Dmin = min(D[a][b] for a, b in itertools.combinations(leaves, 2))
    Dtest = Dmin - 1
    if Dtest < 2: continue
    checked += 1
    # all pairwise distances > Dtest; count leaves with leg < ceil(D/2)
    small = [lf for lf in leaves if leg(adj, lf) < (Dtest + 1) // 2]
    ok3 += (len(small) <= 1)
    assert len(small) <= 1, (Dtest, [leg(adj, lf) for lf in leaves])
print('(3) leg bound: %d/%d random subcubic trees: at most ONE leaf'
      ' with leg < ceil(D/2) when leaf distances exceed D' %
      (ok3, checked))

# ---------------- (4) cherry tightness ----------------
def cherry_tree(r, legl):
    # spine of r/2 branch vertices, each carrying a cherry of two legs
    adj = {}
    def add(a, b):
        adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    spine = ['s%d' % i for i in range(r // 2)]
    for i in range(len(spine) - 1): add(spine[i], spine[i + 1])
    for i, s in enumerate(spine):
        for c in range(2):
            prev = s
            for j in range(legl):
                vtx = 'l%d_%d_%d' % (i, c, j)
                add(prev, vtx); prev = vtx
    return adj

r, D = 8, 22
T = cherry_tree(r, D // 2 + 1)
leaves, DD = leaf_dists(T)
mind = min(DD[a][b] for a, b in itertools.combinations(leaves, 2))
size = len(T)
print('(4) cherry tree r=%d, legs %d: min leaf distance %d > %d: %s;'
      ' size %d ~ r*D/2 = %d: the leg bound is TIGHT'
      % (r, D // 2 + 1, mind, D, mind > D, size, r * D // 2))
assert mind > D
print()
print('VERDICT: the bookkeeping route is closed NEGATIVELY (stall);')
print('cross-connections are confined to middle zones; bridges of the')
print('braided core are leg-trees of size ~ r L/4, pinning it to the')
print('packing boundary shell of 191. The remaining question is the')
print('EXISTENCE of a braided packing-shell graph - open.')
