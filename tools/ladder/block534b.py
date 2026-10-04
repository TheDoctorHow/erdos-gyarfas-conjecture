import itertools
from math import gcd
from collections import deque
from block534 import (matchings, rel, open_legal, girth, ladder_adj,
                      crossing_components, blocks)

# the 34 legal 12-blocks
B12 = [m for m in matchings(list(range(1, 13))) if open_legal(m)]
print("legal 12-blocks:", len(B12))
print("   of these, containing the span-11 chord (1,12):",
      sum(1 for m in B12 if (1, 12) in m))
print()

# ---------- the full disconnected family at 2N = 24 ----------
# a disconnected ladder at 2N=24 is exactly two 12-arcs, each legal.
fam = []
for r in range(12):                       # 12 distinct 12/12 splits
    for m1 in B12:
        for m2 in B12:
            ch = [((a - 1 + r) % 24 + 1, (b - 1 + r) % 24 + 1) for (a, b) in m1]
            ch += [((a + 11 + r) % 24 + 1, (b + 11 + r) % 24 + 1)
                   for (a, b) in m2]
            ch = [(min(e), max(e)) for e in ch]
            fam.append(tuple(sorted(ch)))
fam = sorted(set(fam))
print("candidate disconnected ladders at 2N = 24:", len(fam))

legal = []
for ch in fam:
    adj = ladder_adj(24, ch)
    if adj is None:
        continue
    if girth(adj) >= 5:
        legal.append(ch)
print("   girth >= 5 (BFS-verified, every one):", len(legal))
print("   girth < 5 rejects:", len(fam) - len(legal))

bad = sum(1 for ch in legal if len(crossing_components(list(ch))) < 2)
print("   with connected crossing graph (should be 0):", bad)
print()

# ---------- d distribution: does alignment meet disconnection? ----------
def dval(twoN, ch):
    d = 0
    for (a, b) in ch:
        s = b - a
        d = gcd(gcd(d, s - 1), twoN - s - 1)
    return d

from collections import Counter
dist = Counter(dval(24, ch) for ch in legal)
print("d = gcd(savings) over the disconnected family:")
for k in sorted(dist):
    print("   d = %2d : %d" % (k, dist[k]))
print()

# ---------- the repair: minimal blocks ----------
def minimal_blocks(twoN, ch):
    bl = blocks(twoN, ch)
    out = []
    for (p, L) in bl:
        arc = set(((p - 1 + t) % twoN) + 1 for t in range(L))
        if not any(L2 < L and
                   set(((q - 1 + t) % twoN) + 1
                       for t in range(L2)) <= arc
                   for (q, L2) in bl):
            out.append((p, L))
    return out

def endpoint_lemma_ok(twoN, ch, d):
    """in every minimal block (or the whole ladder if none), the left
    endpoints are not all congruent mod d"""
    mb = minimal_blocks(twoN, ch)
    tests = []
    if not mb:
        tests.append([(min(e), max(e)) for e in ch])
    for (p, L) in mb:
        arc = set(((p - 1 + t) % twoN) + 1 for t in range(L))
        sub = [e for e in ch if e[0] in arc]
        # rotate so the block starts at 1
        sub = [tuple(sorted(((v - p) % twoN) + 1 for v in e)) for e in sub]
        tests.append(sub)
    for sub in tests:
        if len(set(min(e) % d for e in sub)) < 2:
            return False
    return True

def parity_lemma_ok(twoN, ch):
    """in every minimal block with all spans odd, the left endpoints
    take both parities"""
    mb = minimal_blocks(twoN, ch) or [(1, twoN)]
    for (p, L) in mb:
        arc = set(((p - 1 + t) % twoN) + 1 for t in range(L))
        sub = [e for e in ch if e[0] in arc]
        sub = [tuple(sorted(((v - p) % twoN) + 1 for v in e)) for e in sub]
        if all((max(e) - min(e)) % 2 == 1 for e in sub):
            if len(set(min(e) % 2 for e in sub)) < 2:
                return False
    return True

def block_conn_ok(twoN, ch):
    for (p, L) in minimal_blocks(twoN, ch):
        arc = set(((p - 1 + t) % twoN) + 1 for t in range(L))
        sub = [e for e in ch if e[0] in arc]
        if len(crossing_components(sub)) != 1:
            return False
    return True

print("REPAIR CHECKS on the 2N = 24 disconnected family")
print("   minimal blocks have connected crossing graph:",
      all(block_conn_ok(24, ch) for ch in legal), "/", len(legal))
print("   parity lemma holds in every minimal block:",
      all(parity_lemma_ok(24, ch) for ch in legal))
e_ok = all(endpoint_lemma_ok(24, ch, dval(24, ch))
           for ch in legal if dval(24, ch) >= 3)
print("   endpoint lemma holds in every minimal block, d >= 3:", e_ok)
print()

# ---------- the same checks over ALL perfect matchings, 2N <= 12 ----------
print("minimal-block crossing connectivity over ALL matchings")
for n in (8, 10, 12):
    tot = ok = 0
    for m in matchings(list(range(1, n + 1))):
        tot += 1
        if block_conn_ok(n, [(min(e), max(e)) for e in m]):
            ok += 1
    print("   2N = %2d : %6d of %6d" % (n, ok, tot))
