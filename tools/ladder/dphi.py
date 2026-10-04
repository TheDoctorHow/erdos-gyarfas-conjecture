"""The position-sum difference set of the crossing pairs.

phi = x + y.  A crossing pair gives a cycle of length 2 + Dphi (571.2),
so avoiding a cycle of length L means avoiding Dphi = L-2.
"""
import sys, random
from collections import deque
from gwind import lift, girth

def crossing_dphi(mate, n):
    """multiset of Dphi over crossing pairs, x1 taken in [0,n)."""
    isleft, span = lift(mate, n)
    R = [(x, span[x % n]) for x in range(2*n) if isleft[x % n]]
    out = []
    for i in range(len(R)):
        x1, s1 = R[i]
        if x1 >= n: continue
        y1 = x1+s1
        for j in range(i+1, len(R)):
            x2, s2 = R[j]
            if x2 >= y1: break
            y2 = x2+s2
            if y2 > y1: out.append((x2+y2)-(x1+y1))
    return out

def profile(mate, n):
    d = crossing_dphi(mate, n)
    s = set(d)
    return {'girth': girth(mate, n), 'pairs': len(d), 'distinct': len(s),
            'min': min(s) if s else None, 'max': max(s) if s else None,
            'has6': 6 in s, 'has14': 14 in s,
            'mod4': sorted(set(v % 4 for v in s))}

def short(adj, n, a, b, maxlen=4):
    d = {a: 0}; q = deque([a])
    while q:
        u = q.popleft()
        if d[u] >= maxlen-1: continue
        for w in adj[u]:
            if w not in d: d[w] = d[u]+1; q.append(w)
    return b in d

def build_avoiding(n, rng, forbidden, tries=3000):
    """girth-legal closed ladder with no crossing pair whose Dphi is in
    `forbidden`."""
    for _ in range(tries):
        adj = [[(v+1) % n, (v-1) % n] for v in range(n)]
        mate = [-1]*n; ok = True
        placed = []                      # (x, y) with x < y, in [0,n)
        while ok:
            a = next((v for v in range(n) if mate[v] == -1), None)
            if a is None: break
            cand = []
            for b in range(n):
                if b == a or mate[b] != -1: continue
                if short(adj, n, a, b): continue
                trial = placed + [(min(a, b), max(a, b))]
                if not crossing_bad(trial, n, forbidden): cand.append(b)
            if not cand: ok = False; break
            b = rng.choice(cand)
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
            placed.append((min(a, b), max(a, b)))
        if ok and all(m != -1 for m in mate):
            return mate
    return None

def crossing_bad(chords, n, forbidden):
    """any crossing pair among these chords (each given as residues) with
    Dphi in forbidden?  Lift each chord to its short form."""
    R = []
    for (u, v) in chords:
        d = v-u
        if d <= n-d: R.append((u, d))
        else: R.append((v, n-d))
    R2 = R + [(x+n, s) for x, s in R]
    R2.sort()
    for i in range(len(R2)):
        x1, s1 = R2[i]
        if x1 >= n: continue
        y1 = x1+s1
        for j in range(i+1, len(R2)):
            x2, s2 = R2[j]
            if x2 >= y1: break
            y2 = x2+s2
            if y2 > y1 and ((x2+y2)-(x1+y1)) in forbidden: return True
    return False

if __name__ == '__main__':
    from horizon import try_build
    rng = random.Random(24680)
    print("== octagon-free ladders: their Dphi structure ==")
    for n in (28, 34, 38, 44):
        m = try_build(n, rng, forbid=(8,), tries=3000)
        if m is None: print("  2N=%d none" % n); continue
        p = profile(m, n)
        print("  2N=%2d  girth %d  crossing pairs %4d  Dphi in [%d,%d]  "
              "6? %-5s 14? %-5s  residues mod 4: %s"
              % (n, p['girth'], p['pairs'], p['min'], p['max'],
                 p['has6'], p['has14'], p['mod4']))
        sys.stdout.flush()
