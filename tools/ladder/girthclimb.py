"""Maximise the girth of a closed ladder by 2-swaps on the matching."""
import sys, random, time
from collections import deque

def girth_and_witness(mate, n, cap=10**9):
    A = [[(v+1) % n, (v-1) % n, mate[v]] for v in range(n)]
    best = cap; wit = None
    for src in range(n):
        d = {src: 0}; par = {src: -1}; q = deque([src])
        while q:
            u = q.popleft()
            if 2*d[u] >= best: break
            for w in A[u]:
                if w not in d:
                    d[w] = d[u]+1; par[w] = u; q.append(w)
                elif w != par[u]:
                    L = d[u]+d[w]+1
                    if L < best:
                        best = L
                        # recover the cycle
                        pu = []; x = u
                        while x != -1: pu.append(x); x = par[x]
                        pw = []; x = w
                        while x != -1: pw.append(x); x = par[x]
                        wit = pu + pw[::-1]
    return best, wit

def swap(mate, a, c):
    b = mate[a]; d = mate[c]
    if len({a, b, c, d}) < 4: return False
    mate[a] = c; mate[c] = a; mate[b] = d; mate[d] = b
    return True

def climb(n, rng, seconds=60.0, restarts=6):
    bestg = 0; bestm = None
    t0 = time.time()
    for _ in range(restarts):
        # random perfect matching on the cycle
        pts = list(range(n)); rng.shuffle(pts)
        mate = [-1]*n
        for i in range(0, n, 2):
            a, b = pts[i], pts[i+1]
            mate[a] = b; mate[b] = a
        g, wit = girth_and_witness(mate, n)
        while time.time()-t0 < seconds:
            if wit is None: break
            # pick a chord on the shortest cycle, swap it with a random chord
            chords = [v for v in wit if mate[v] in wit]
            a = rng.choice(chords if chords else wit)
            c = rng.randrange(n)
            keep = (mate[a], mate[c])
            if not swap(mate, a, c): continue
            g2, wit2 = girth_and_witness(mate, n)
            if g2 >= g:
                g, wit = g2, wit2
            else:
                swap(mate, a, keep[0])       # undo
                mate[a] = keep[0]; mate[keep[0]] = a
                mate[c] = keep[1]; mate[keep[1]] = c
            if g > bestg:
                bestg = g; bestm = list(mate)
        if time.time()-t0 >= seconds: break
    return bestg, bestm

if __name__ == '__main__':
    rng = random.Random(777)
    print("2N    best girth found")
    for n in (40, 60, 80, 100, 140, 200, 300):
        g, m = climb(n, rng, seconds=45.0, restarts=4)
        print("  %4d      %2d" % (n, g))
        if m is not None:
            import json, io
            io.open('ladder_%d.json' % n, 'w').write(json.dumps(m))
        sys.stdout.flush()
