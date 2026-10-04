"""Block every k <= 2 mechanism for BOTH 8 and 16 at once:
   no span in {7,15}, no nested pair with Dsigma in {6,14},
   no crossing pair with Dphi in {6,14}.
Which ladders survive, and do they still have an 8- or 16-cycle?
"""
import sys, random
from collections import deque
from gwind import lift, girth

TARGETS = (8, 16)
SPANS = set(L-1 for L in TARGETS)
DIFFS = set(L-2 for L in TARGETS)

def short(adj, n, a, b, maxlen=4):
    d = {a: 0}; q = deque([a])
    while q:
        u = q.popleft()
        if d[u] >= maxlen-1: continue
        for w in adj[u]:
            if w not in d: d[w] = d[u]+1; q.append(w)
    return b in d

def to_rungs(chords, n):
    R = []
    for (u, v) in chords:
        d = v-u
        if d <= n-d: R.append((u, d))
        else: R.append((v, n-d))
    return R

def blocked(chords, n):
    """True if some k<=2 mechanism fires for 8 or 16."""
    R = to_rungs(chords, n)
    for x, s in R:
        if s in SPANS: return True
    R2 = sorted(R + [(x+n, s) for x, s in R])
    for i in range(len(R2)):
        x1, s1 = R2[i]
        if x1 >= n: continue
        y1 = x1+s1
        for j in range(i+1, len(R2)):
            x2, s2 = R2[j]
            if x2 >= y1: break
            y2 = x2+s2
            if y2 > y1:
                if ((x2+y2)-(x1+y1)) in DIFFS: return True
            elif y2 < y1:
                if (s1-s2) in DIFFS: return True
    return False

def build(n, rng, tries=4000):
    for _ in range(tries):
        adj = [[(v+1) % n, (v-1) % n] for v in range(n)]
        mate = [-1]*n; placed = []; ok = True
        while ok:
            a = next((v for v in range(n) if mate[v] == -1), None)
            if a is None: break
            cand = []
            for b in range(n):
                if b == a or mate[b] != -1: continue
                if short(adj, n, a, b): continue
                if not blocked(placed + [(min(a, b), max(a, b))], n):
                    cand.append(b)
            if not cand: ok = False; break
            b = rng.choice(cand)
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
            placed.append((min(a, b), max(a, b)))
        if ok and all(m != -1 for m in mate):
            return mate
    return None

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    from no8 import has_len
    rng = random.Random(86420)
    print("ladders with EVERY k<=2 mechanism blocked for both 8 and 16:")
    tot = 0; still = 0
    for n in range(14, 61, 2):
        m = build(n, rng, tries=600)
        if m is None:
            print("  2N=%2d : none in 600 attempts" % n); sys.stdout.flush(); continue
        tot += 1
        isleft, span = lift(m, n)
        sp = sorted(set(span[r] for r in range(n) if isleft[r]))
        h8 = has_len(m, n, 8); h16 = has_len(m, n, 16)
        if h8 or h16: still += 1
        print("  2N=%2d : FOUND  girth %d  spans %s   8-cycle %-5s 16-cycle %-5s"
              % (n, girth(m, n), sp if len(sp) < 8 else '%d values' % len(sp),
                 h8, h16))
        sys.stdout.flush()
    print()
    print("built %d ; still have an 8- or 16-cycle (so from k >= 3): %d"
          % (tot, still))
