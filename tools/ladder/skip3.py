import sys, random
from collections import deque

# ---------------- ladders ----------------

def short_cycle(mate, n, a, b, maxlen=4):
    A = [[(v+1) % n, (v-1) % n] for v in range(n)]
    for v in range(n):
        if mate[v] is not None: A[v].append(mate[v])
    dist = {a: 0}; q = deque([a])
    while q:
        u = q.popleft()
        if dist[u] >= maxlen-1: continue
        for w in A[u]:
            if w not in dist:
                dist[w] = dist[u]+1; q.append(w)
    return b in dist

def all_ladders(n):
    mate = [None]*n; out = []
    def rec():
        a = next((v for v in range(n) if mate[v] is None), None)
        if a is None:
            out.append(list(mate)); return
        for b in range(a+1, n):
            if mate[b] is None and not short_cycle(mate, n, a, b):
                mate[a] = b; mate[b] = a
                rec()
                mate[a] = None; mate[b] = None
    rec()
    return out

def random_ladder(n, rng, tries=4000):
    for _ in range(tries):
        mate = [None]*n; ok = True
        while ok:
            a = next((v for v in range(n) if mate[v] is None), None)
            if a is None: break
            cand = [b for b in range(a+1, n)
                    if mate[b] is None and not short_cycle(mate, n, a, b)]
            if not cand: ok = False; break
            b = rng.choice(cand)
            mate[a] = b; mate[b] = a
        if ok and all(m is not None for m in mate):
            return mate
    return None

# ------------- the one-lap family, by DFS -------------

def onelap_lengths(mate, n):
    """set of lengths of winding-1 all-forward simple cycles."""
    d = [(mate[v]-v) % n for v in range(n)]
    lens = set()

    def bb(p, c):                      # overshoot of chord p past c
        return max(0, d[p] - ((c-p) % n))

    for x1 in range(n):
        # k = 1
        lens.add(n - d[x1] + 1)
        for x2 in range(x1+1, n):
            g = x2 - x1
            if g == d[x1]: continue
            b1 = max(0, d[x1]-g)
            # k = 2 : close straight back
            g2 = n - g
            if g2 != d[x2]:
                b2 = max(0, d[x2]-g2)
                if b1+b2 < d[x2] and b2+b1 < d[x1]:
                    lens.add(n - (d[x1]-1) - (d[x2]-1) + 2*(b1+b2))
            # k >= 3
            t1 = d[x1]-1-2*b1
            stack = [(x1, x2, b1, t1)]
            while stack:
                p, c, bp, tot = stack.pop()
                dc = d[c]
                for nx in range(c+1, n):
                    g = nx - c
                    if g == dc: continue
                    bc = max(0, dc-g)
                    if bp + bc >= dc: continue
                    stack.append((c, nx, bc, tot + dc-1-2*bc))
                # close: next is x1
                gcl = n - c + x1
                if gcl != dc:
                    bcl = max(0, dc-gcl)
                    if bp + bcl < dc and bcl + b1 < d[x1]:
                        lens.add(n - (tot + dc-1-2*bcl))
    return lens

# ------------- ground truth -------------

def adj_of(mate, n):
    return [[(v+1) % n, (v-1) % n, mate[v]] for v in range(n)]

def all_cycle_lengths(mate, n):
    A = adj_of(mate, n); lens = set()
    def rec(start, u, vis, ln):
        for w in A[u]:
            if w == start and ln >= 3: lens.add(ln)
            elif w > start and not (vis >> w) & 1:
                rec(start, w, vis | (1 << w), ln+1)
    for s0 in range(n):
        rec(s0, s0, 1 << s0, 1)
    return lens

def girth(mate, n):
    A = adj_of(mate, n); g = 10**9
    for src in range(n):
        dist = {src: 0}; par = {src: -1}; q = deque([src])
        while q:
            u = q.popleft()
            for w in A[u]:
                if w not in dist:
                    dist[w] = dist[u]+1; par[w] = u; q.append(w)
                elif w != par[u]: g = min(g, dist[u]+dist[w]+1)
    return g

def pow2(S):
    return sorted(x for x in S if x >= 4 and (x & (x-1)) == 0)

# ------------- driver -------------

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'exact':
        for n in [12, 14, 16, 18]:
            Ls = [m for m in all_ladders(n) if girth(m, n) >= 5]
            miss = 0; unsound = 0; ex = []
            for mate in Ls:
                one = onelap_lengths(mate, n)
                full = all_cycle_lengths(mate, n)
                if not (one <= full):
                    unsound += 1
                    if unsound == 1: print("  UNSOUND", mate, sorted(one-full))
                if not pow2(one):
                    miss += 1
                    if len(ex) < 3: ex.append((list(mate), sorted(full), sorted(one)))
            print("2N=%2d  ladders %6d   one-lap misses a power of two: %d"
                  "   unsound: %d" % (n, len(Ls), miss, unsound))
            for e in ex:
                print("    mate %s" % (e[0],))
                print("      full    %s" % (e[1],))
                print("      one-lap %s" % (e[2],))
            sys.stdout.flush()
    else:
        rng = random.Random(20260810)
        for n in [20, 22, 24, 26, 28, 30]:
            miss = 0; cnt = 0; ex = []
            for _ in range(300):
                mate = random_ladder(n, rng)
                if mate is None or girth(mate, n) < 5: continue
                cnt += 1
                one = onelap_lengths(mate, n)
                if not pow2(one):
                    miss += 1
                    if len(ex) < 2: ex.append((list(mate), sorted(one)))
            print("2N=%2d  sampled %4d   one-lap misses a power of two: %d"
                  % (n, cnt, miss))
            for e in ex:
                print("    mate %s" % (e[0],))
                print("      one-lap %s" % (e[1],))
            sys.stdout.flush()
