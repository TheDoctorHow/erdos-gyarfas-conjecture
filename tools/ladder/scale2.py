"""What girth forces on the chord diagram, and how far three chords reach."""
import sys, random
from collections import deque
from gwind import lift, girth
from scale import triple_lengths, pows

# ---------- the girth constraints, checked ----------

def check_constraints(mate, n):
    """single: sigma >= g-1 ; nested: sigma_out - sigma_in >= g-2 ;
    crossing: 2*dx + (sigma2 - sigma1) >= g-2."""
    g = girth(mate, n)
    isleft, span = lift(mate, n)
    R = [(x, span[x % n]) for x in range(2*n) if isleft[x % n]]
    bad = []
    for x, s in R:
        if x < n and s < g-1: bad.append(('single', x, s))
    for i in range(len(R)):
        x1, s1 = R[i]
        if x1 >= n: continue
        y1 = x1+s1
        for j in range(i+1, len(R)):
            x2, s2 = R[j]
            if x2 >= y1: break
            y2 = x2+s2
            if y2 < y1:                       # nested, 1 outer
                if s1-s2 < g-2: bad.append(('nested', x1, s1, x2, s2))
            elif y2 > y1:                     # crossing
                if 2*(x2-x1)+(s2-s1) < g-2:
                    bad.append(('crossing', x1, s1, x2, s2))
    return g, bad

# ---------- structured ladders ----------

def short_cycle(mate, n, a, b, maxlen=4):
    A = [[(v+1) % n, (v-1) % n] for v in range(n)]
    for v in range(n):
        if mate[v] is not None: A[v].append(mate[v])
    d = {a: 0}; q = deque([a])
    while q:
        u = q.popleft()
        if d[u] >= maxlen-1: continue
        for w in A[u]:
            if w not in d: d[w] = d[u]+1; q.append(w)
    return b in d

def ladder_with_spans(n, rng, allowed, tries=6000):
    """girth-legal closed ladder whose circular spans all lie in `allowed`."""
    for _ in range(tries):
        mate = [None]*n; ok = True
        while ok:
            a = next((v for v in range(n) if mate[v] is None), None)
            if a is None: break
            cand = [b for b in range(a+1, n)
                    if mate[b] is None and min((b-a) % n, (a-b) % n) in allowed
                    and not short_cycle(mate, n, a, b)]
            if not cand: ok = False; break
            b = rng.choice(cand); mate[a] = b; mate[b] = a
        if ok and all(m is not None for m in mate) and girth(mate, n) >= 5:
            return mate
    return None

if __name__ == '__main__':
    from highgirth import build
    from skip3 import random_ladder
    rng = random.Random(31337)

    print("== the girth constraints on the chord diagram ==")
    tested = 0; viol = 0
    for n in (16, 20, 24, 30, 40, 60):
        for _ in range(12):
            m = random_ladder(n, rng)
            if m is None or girth(m, n) < 5: continue
            g, bad = check_constraints(m, n)
            tested += 1; viol += len(bad)
            if bad: print("  VIOLATION 2N=%d g=%d %s" % (n, g, bad[:2]))
    for n, tg in ((40, 7), (80, 8)):
        m = build(n, rng, tg, tries=40)
        if m:
            g, bad = check_constraints(m, n)
            tested += 1; viol += len(bad)
            if bad: print("  VIOLATION high-girth 2N=%d g=%d %s" % (n, g, bad[:2]))
    print("  ladders checked %d   violations %d" % (tested, viol))

    print()
    print("== how far do three chords reach?  (span diversity is the fuel) ==")
    def row(m, n, tag):
        g = girth(m, n)
        isleft, span = lift(m, n)
        sp = sorted(set(span[r] for r in range(n) if isleft[r]))
        T = triple_lengths(m, n)
        ks = sorted(T)
        hit = [L for L in pows(n) if L in T]
        print("  %-26s 2N=%3d g=%d  distinct spans %2d  triples give %3d "
              "lengths, top %3d (%.2f of 2N)  pow2 %s"
              % (tag, n, g, len(sp), len(ks), ks[-1] if ks else -1,
                 (ks[-1]/n) if ks else 0, hit))
        return bool(hit)
    fails = 0; tot = 0
    for n in (24, 32, 40, 60):
        for lbl, allowed in (("2 spans", None), ("3 spans", None)):
            k = 2 if lbl == "2 spans" else 3
            allowed = set(rng.sample(range(4, n//2+1), k))
            m = ladder_with_spans(n, rng, allowed)
            if m is None:
                print("  %-26s 2N=%3d  none exists with spans %s"
                      % (lbl, n, sorted(allowed)))
                continue
            tot += 1
            if not row(m, n, "%s %s" % (lbl, sorted(allowed))): fails += 1
        sys.stdout.flush()
    print("  structured ladders tested %d   three chords failed on %d"
          % (tot, fails))
