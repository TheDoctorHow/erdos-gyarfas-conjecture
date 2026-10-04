"""Exhaustive DFS over attachment matchings of the l=11 K3,3 template.

SOUNDNESS: every filter below is a PROVEN-NECESSARY condition for a
mod-8-free girth->=7 completion (one-edge tables + 2-web tables derived
in PROOF-PROGRAM.md section 5 chain). Filters only ever REJECT pairings
that provably create a cycle of length ==0 (mod 8) or girth < 7.
If the search exhausts with no survivor, no completion exists at l=11.
Survivors (if any) are written out for full oracle verification.
"""
import sys
sys.setrecursionlimit(10000)

L = 11  # corridor length, offsets 1..10

# corridor id c = 0..8 encodes (x, y) = (c//3, c%3); u-end u_x, v-end v_y.
def rel(c1, c2):
    """'same' / 'adjU' / 'adjV' / 'disj' for two corridors."""
    if c1 == c2:
        return 'same'
    x1, y1 = divmod(c1, 3)
    x2, y2 = divmod(c2, 3)
    if x1 == x2:
        return 'adjU'
    if y1 == y2:
        return 'adjV'
    return 'disj'

def edge_ok(c1, k1, c2, k2):
    """One-edge (proven-necessary) tables. k = u-offset 1..10."""
    r = rel(c1, c2)
    if r == 'same':
        d = abs(k1 - k2)
        return d >= 6 and d % 8 not in (5, 7)
    if r == 'adjU':
        i, j = k1, k2                    # offsets from shared u-vertex
    elif r == 'adjV':
        i, j = L - k1, L - k2            # offsets from shared v-vertex
    else:
        # disjoint: parity law + route cycles (canonical u-frames)
        a, b = k1, k2
        if (a - b) % 2:
            return False
        if (a - b) % 8 in (1, 7):        # 23 +/- (a-b) == 0 (mod 8)
            return False
        if (a + b) % 8 in (1, 5):        # 23+(a+b), 45-(a+b)
            return False
        return True
    s = i + j
    if s < 6:
        return False                      # through-cycle girth
    if s % 8 in (3, 5, 7):
        return False
    if (i - j) % 8 in (3, 5) or (j - i) % 8 in (3, 5):
        return False
    return True

def pair_ok(e, f):
    """2-web (proven-necessary) tables for two placed attachment edges."""
    (c1, k1, c2, k2) = e
    (c3, k3, c4, k4) = f
    # same-corridor internal pairs
    if c1 == c2 and c3 == c4 and c1 == c3:
        x1, y1 = min(k1, k2), max(k1, k2)
        x2, y2 = min(k3, k4), max(k3, k4)
        if x2 < x1:
            x1, y1, x2, y2 = x2, y2, x1, y1
        if x2 > y1:                       # disjoint spans
            if ((y1 - x1) + (y2 - x2)) % 8 == 6:
                return False
        elif y2 < y1:                     # nested
            s = (x2 - x1) + (y1 - y2)
            if s % 8 == 6 or s < 5:
                return False
        else:                             # crossing
            s = (x2 - x1) + (y2 - y1)
            if s % 8 == 6 or s < 5:
                return False
        return True
    # internal edge on c + internal edge on c' (2-around web)
    if c1 == c2 and c3 == c4:
        d1 = abs(k1 - k2)
        d2 = abs(k3 - k4)
        r = rel(c1, c3)
        if r in ('adjU', 'adjV'):
            if (d1 + d2) % 8 in (4, 6):
                return False
        else:
            if (d1 + d2) % 8 == 4:
                return False
        return True
    # same-pair quad: both edges between the same corridor pair {A, B}:
    # near-near cycle 2 + |dk_A| + |dk_B| must have length >= 7 and != 0 mod 8
    pa, pb = frozenset((c1, c2)), frozenset((c3, c4))
    if pa == pb and len(pa) == 2:
        A = min(c1, c2)
        a1 = k1 if c1 == A else k2
        b1 = k2 if c1 == A else k1
        a2 = k3 if c3 == A else k4
        b2 = k4 if c3 == A else k3
        s = abs(a1 - a2) + abs(b1 - b2)
        if s <= 4 or s % 8 == 6:
            return False
    return True   # remaining cross-geometry tables omitted (necessary-only)

def feet(e):
    (c1, k1, c2, k2) = e
    return [(c1, k1), (c2, k2)]

def search():
    order = [(c, k) for c in range(9) for k in range(1, L)]
    matched = {}
    edges = []
    stats = {'nodes': 0, 'sols': 0}
    cap = 2 * 10**9

    def bt(idx):
        stats['nodes'] += 1
        if stats['nodes'] % 50_000_000 == 0:
            print('nodes=%dM depth=%d edges=%d' %
                  (stats['nodes'] // 10**6, idx, len(edges)), flush=True)
        if stats['nodes'] > cap:
            raise SystemExit('NODE CAP HIT: INCONCLUSIVE')
        while idx < len(order) and order[idx] in matched:
            idx += 1
        if idx == len(order):
            stats['sols'] += 1
            print('SURVIVOR:', edges, flush=True)
            with open('dfs11-survivors.txt', 'a') as f:
                f.write(repr(edges) + '\n')
            return
        c1, k1 = order[idx]
        for j in range(idx + 1, len(order)):
            c2, k2 = order[j]
            if (c2, k2) in matched:
                continue
            if not edge_ok(c1, k1, c2, k2):
                continue
            e = (c1, k1, c2, k2)
            if all(pair_ok(e, f) for f in edges):
                matched[(c1, k1)] = matched[(c2, k2)] = True
                edges.append(e)
                bt(idx + 1)
                edges.pop()
                del matched[(c1, k1)], matched[(c2, k2)]
        return

    bt(0)
    print('EXHAUSTED: nodes=%d survivors=%d' % (stats['nodes'], stats['sols']),
          flush=True)

if __name__ == '__main__':
    search()
