import itertools
from math import gcd
from block534 import (matchings, rel, open_legal, girth, ladder_adj,
                      crossing_components, blocks, spectrum)

def aligned(m, d):
    return all((b - a) % d == 1 for (a, b) in m)

B12 = [m for m in matchings(list(range(1, 13)))
       if open_legal(m) and aligned(m, 3)]
B14 = [m for m in matchings(list(range(1, 15)))
       if open_legal(m) and aligned(m, 3)]
print("legal blocks with every span = 1 mod 3:  size 12:", len(B12),
      "  size 14:", len(B14))
for m in B12:
    print("   12-block", m, sorted(b - a for a, b in m))
for m in B14[:4]:
    print("   14-block", m, sorted(b - a for a, b in m))
print()

# ---------- assemble 2N = 26 : block of 12 then block of 14 ----------
best = None
for m1 in B12:
    for m2 in B14:
        ch = sorted((min(e), max(e)) for e in
                    m1 + [(a + 12, b + 12) for (a, b) in m2])
        adj = ladder_adj(26, ch)
        if adj is None or girth(adj) < 5:
            continue
        d = 0
        for (a, b) in ch:
            s = b - a
            d = gcd(gcd(d, s - 1), 26 - s - 1)
        if d >= 3 and len(crossing_components(ch)) > 1:
            best = (ch, d, adj)
            break
    if best:
        break

if not best:
    print("no disconnected aligned ladder at 2N = 26")
else:
    ch, d, adj = best
    print("DISCONNECTED **AND** ALIGNED at 2N = 26")
    print("   chords:", ch)
    print("   girth (BFS) =", girth(adj), "  degrees =",
          sorted(set(len(a) for a in adj)))
    print("   d = gcd(savings) =", d)
    comps = crossing_components(ch)
    print("   crossing components:", len(comps))
    for c in comps:
        print("      ", [ch[i] for i in c])
    print("   proper self-matched arcs:", blocks(26, ch))

    # the shortcut mechanism really is confined
    print()
    print("   2N mod d = 26 mod %d = %d" % (d, 26 % d))
    sp = spectrum(adj)
    print("   full cycle spectrum:", sp)
    print("   residues mod %d realised:" % d, sorted(set(x % d for x in sp)))
    print("   powers of two present:", [x for x in sp if x & (x - 1) == 0])

    # what the OLD theorem could say: nothing, hypothesis fails.
    # what the REPAIRED theorem says: work inside a minimal block.
    for (p, L) in blocks(26, ch):
        arc = set(((p - 1 + t) % 26) + 1 for t in range(L))
        sub = [e for e in ch if e[0] in arc]
        sub_r = [tuple(sorted(((v - p) % 26) + 1 for v in e)) for e in sub]
        lefts = sorted(min(e) % d for e in sub_r)
        print("   block (%d,%d): own components %d, left endpoints mod %d = %s"
              % (p, L, len(crossing_components(sub)), d, lefts))
        # the crossing pair the repaired proof produces
        for i, j in itertools.combinations(range(len(sub_r)), 2):
            if rel(sub_r[i], sub_r[j]) == 'X':
                da = min(sub_r[j]) - min(sub_r[i])
                if da % d:
                    L1 = abs(min(sub_r[j]) - min(sub_r[i])) + \
                         abs(max(sub_r[j]) - max(sub_r[i]))
                    print("      escaping crossing pair", sub[i], sub[j],
                          " cycle lengths", 2 + L1, "and", 26 + 2 - L1,
                          " residues", (2 + L1) % d, (28 - L1) % d,
                          " vs shortcut class", 26 % d)
                    break
