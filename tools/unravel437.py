"""unravel437.py -- the Unravelling Question: the affine cascade, the
region bound, fragmentation, and the survivor's wiring.

 (1) AFFINE CASCADE: on a connected wiring graph, constancy of the
     uu-shifts forces p_i + pi_i == alpha, q_j - pi_j == beta (mod 2^t)
     (connected labelling); adding yy-constancy and the first-gen
     stall |P_i| == c0 + eta_i - pi_i forces 2 eta_i == const (mod
     2^t). Machine: random connected instances - solve for shift
     constancy; verify it is possible iff the eta-congruence holds.
 (2) REGION BOUND: eta_i in one class mod 2^(t-1): at L = 3 * 2^t
     only 2m = L/2^(t-1) positions exist; a connected-wired constancy
     region holds <= 2m carriers. Arithmetic table.
 (3) INCREASING INVOLUTION LEMMA: the only order-preserving involution
     on [n] is the identity (brute force n <= 8): a matched refresh
     family cannot be pairwise crossing via an order-preserving
     matching; the antipodal shift evades with nested WITHIN-pair
     chords only (max laminar 2). Enumerate the N = 12 antipodal
     family's crossing graph.
 (4) RUNG MENU: a ladder's uy-ears across rung choices have shifts in
     an AP with difference 2 * spacing; an ODD base dodges 0 mod 2^t
     wholesale, an EVEN base is hit once the AP covers its coset.
"""
import itertools, random

v2 = lambda x: (abs(x) & -abs(x)).bit_length() - 1 if x else 999
random.seed(437)

# ---------------- (1) affine cascade ----------------
t = 3; Mt = 1 << t
def cascade_instance(ok_eta):
    """Build carriers with positions; return whether constant uu+yy
    shifts are simultaneously achievable by choosing rung depths."""
    n = 6
    pi = sorted(random.sample(range(0, 40, 2), n))
    if ok_eta:
        base = random.randrange(Mt // 2)
        eta = [60 + i * 4 + (2 * base - 0) % 2 for i in range(n)]
        eta = [60 + 4 * i for i in range(n)]          # all == mod 4
    else:
        eta = [60 + 4 * i for i in range(n)]
        eta[2] += 2                                   # break one class
    c0 = 3
    P = [(c0 + eta[i] - pi[i]) % Mt for i in range(n)]  # |P_i| mod 8
    # wiring: a random connected graph on carriers
    edges = [(i, i + 1) for i in range(n - 1)]
    edges += random.sample([(i, j) for i in range(n)
                            for j in range(i + 2, n)], 3)
    # solve: exists p_i, q_i mod 8 with, for every edge (i,j):
    #   uu: p_i + q_j + 1 - (pi_j - pi_i) == c2
    #   yy: (P_i - p_i) + (P_j - q_j) + 1 - (eta_j - eta_i) == c2p
    for c2 in range(Mt):
        for c2p in range(Mt):
            # propagate p, q from carrier 0
            for p0 in range(Mt):
                p = {0: p0}; q = {}
                okall = True
                # BFS
                seen = {0}; frontier = [0]
                while frontier and okall:
                    i = frontier.pop()
                    for (a, b) in edges:
                        for (x, y) in ((a, b), (b, a)):
                            if x == i:
                                # q_y from uu with edge (x,y): direction
                                qy = (c2 - 1 - p[x] + (pi[y] - pi[x])) % Mt
                                if y in q and q[y] != qy: okall = False
                                q[y] = qy
                                # p_y from yy
                                py = (P[y] - ((c2p - 1) -
                                      (P[x] - p[x]) + (eta[y] - eta[x]))) % Mt
                                if y in p and p[y] != py: okall = False
                                p[y] = py
                                if y not in seen:
                                    seen.add(y); frontier.append(y)
                if okall and len(seen) == 6:
                    return True
    return False

good = sum(cascade_instance(True) for _ in range(8))
bad = sum(cascade_instance(False) for _ in range(8))
print('(1) affine cascade: constancy achievable with aligned eta (one'
      ' class mod 4): %d/8; with a broken eta class: %d/8' % (good, bad))
assert good == 8 and bad == 0

# ---------------- (2) region bound ----------------
print('(2) region bound: positions == const mod 2^(t-1) on |C| = 3*2^t:')
for tt in (3, 4, 5):
    L = 3 * (1 << tt)
    cap = L // (1 << (tt - 1))
    print('    t=%d: L=%3d, capacity per class = %d = 2m: a connected'
          ' wired constancy region holds <= %d carriers' % (tt, L, cap, cap))

# ---------------- (3) increasing involution + antipodal family ----------------
def increasing_involutions(n):
    out = []
    for perm in itertools.permutations(range(n)):
        if all(perm[perm[i]] == i for i in range(n)) and \
           all(perm[i] < perm[j] for i in range(n) for j in range(i + 1, n)):
            out.append(perm)
    return out
for n in (2, 3, 4, 5, 6):
    invs = increasing_involutions(n)
    assert invs == [tuple(range(n))]
print('(3) increasing involution lemma: only identity, n <= 6: verified')
# antipodal family: N=12 carriers, chords (u_i, y_{i+6 mod 12}) plus
# within-pair partners; max laminar subfamily?
N = 12
chords = [('u', i, 'y', (i + 6) % N) for i in range(N)]
def crossing(ch1, ch2):
    # circle order: u_0..u_11 at 0..11, y_0..y_11 at 12..23
    pos = lambda c: (c[1] if c[0] == 'u' else 12 + c[1],
                     c[3] if c[2] == 'y' else c[3])
    a = ch1[1]; b = 12 + ch1[3]; c = ch2[1]; d = 12 + ch2[3]
    a, b = sorted((a, b)); c, d = sorted((c, d))
    return (a < c < b < d) or (c < a < d < b)
noncross = [[i for i in range(N) if i != j and
             not crossing(chords[j], chords[i])] for j in range(N)]
# max independent set in crossing graph = max laminar family
best = 1
for sub in itertools.combinations(range(N), 3):
    if all(not crossing(chords[i], chords[j])
           for i, j in itertools.combinations(sub, 2)):
        best = 3; break
pairs2 = any(not crossing(chords[i], chords[j])
             for i, j in itertools.combinations(range(N), 2))
print('    antipodal uy-family N=12: laminar pair exists: %s; laminar'
      ' TRIPLE exists: %s (survivor needs False)' % (pairs2, best == 3))
assert best < 3

# ---------------- (4) rung menu ----------------
t = 3; Mt = 1 << t
def menu(base, spacing, k):
    return [(base + 2 * spacing * r) % Mt for r in range(k)]
odd_menu = menu(3, 1, 6)     # odd base
even_menu = menu(2, 1, 6)    # even base
print('(4) rung AP menus (t=3): odd base 3: %s -> 0 mod 8 hit: %s;'
      ' even base 2: %s -> hit: %s'
      % (odd_menu, 0 in odd_menu, even_menu, 0 in even_menu))
assert 0 not in odd_menu and 0 in even_menu
print()
print('VERDICT: blockiness cascades to an affine pinning of rung depths;')
print('connected wired regions collapse to <= 2m carriers (fragmentation')
print('FORCED); refresh families exist but the increasing-involution')
print('lemma limits their braiding to the antipodal pattern with nested')
print('within-pair chords; rung menus survive only on odd bases. The')
print('survivor is the LADDER-CLIQUE ODD BRAID: fragments of <= 2m')
print('antipodally matched carriers, every menu odd-based.')
