import sys, itertools, random

def uniform(n, s):
    # vertices 0..n-1 correspond to 1..n ; chords {i,i+s} for i odd (1-based)
    mate = [None]*n
    for i in range(1, n+1, 2):
        a = i-1
        b = (i+s-1) % n
        mate[a] = b
        mate[b] = a
    assert all(m is not None for m in mate)
    return mate

def walk(mate, n, S):
    """S = sorted list of chord tails.  Canonical one-lap walk.
    Returns (seq, ok) where ok False if some run has length 0."""
    k = len(S)
    seq = []
    for j in range(k):
        x = S[j]; y = S[(j+1) % k]
        d = (mate[x]-x) % n
        gap = (y-x) % n
        if k == 1: gap = n
        if gap == d:
            return None
        h = (x+d) % n
        seq.append(x); seq.append(h)
        r = gap - d
        if r > 0:
            for t in range(1, r): seq.append((h+t) % n)
        else:
            for t in range(1, -r): seq.append((h-t) % n)
    return seq

def criterion(mate, n, S):
    k = len(S)
    d = [ (mate[S[j]]-S[j]) % n for j in range(k) ]
    gap = [ (S[(j+1)%k]-S[j]) % n for j in range(k) ]
    if k == 1: gap = [n]
    if any(g == dd for g, dd in zip(gap, d)): return False
    b = [ max(0, d[j]-gap[j]) for j in range(k) ]
    for j in range(k):
        if b[j] + b[(j+1)%k] >= d[(j+1)%k]:
            return False
    return True

def issimple(seq):
    return seq is not None and len(set(seq)) == len(seq)

# ---------- 1. the 564.0 counterexample ----------
print("== 1. counterexample to 564.0 as stated ==")
n, s = 24, 5
mate = uniform(n, s)
# 1-based starts 1,3, then forward runs of 5: walk given by runs -3,-3,+5,+5
# reconstruct the tails
x = 0  # vertex 1
tails = [x]
runs = [(-1,3), (-1,3), (1,5), (1,5)]
for eps, r in runs[:-1]:
    x = (x + s + eps*r) % n
    tails.append(x)
print("  tails (1-based):", [t+1 for t in sorted(tails)])
S = sorted(tails)
seq = walk(mate, n, S)
print("  displacement check  sum(s+eps*r) =", sum(s+e*r for e, r in runs), " 2N =", n)
print("  walk (1-based):", [v+1 for v in seq])
print("  distinct:", len(set(seq)), " length:", len(seq), " SIMPLE:", issimple(seq))
print("  criterion says:", criterion(mate, n, S))

# ---------- 2. exhaustive: criterion == simplicity, uniform ladders ----------
print()
print("== 2. exhaustive check on uniform ladders ==")
bad = 0; tot = 0; agree = 0
for n in range(12, 27, 2):
    N = n//2
    for s in range(5, N-1, 2):
        mate = uniform(n, s)
        verts = list(range(n))
        for k in range(1, 7):
            for S in itertools.combinations(verts, k):
                S = list(S)
                seq = walk(mate, n, S)
                c = criterion(mate, n, S)
                sm = issimple(seq)
                tot += 1
                if c == sm: agree += 1
                else:
                    bad += 1
                    if bad <= 5:
                        print("   MISMATCH n=%d s=%d S=%s crit=%s simple=%s"
                              % (n, s, S, c, sm))
print("  subsets tested %d   agreements %d   mismatches %d" % (tot, agree, bad))
