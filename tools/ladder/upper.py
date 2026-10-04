"""The matching upper bound for uniform ladders.

Deficiency generators in U(2N,s):
   one disjoint chord          s - 1
   one crossing pair, gap Dx   2 Dx - 2,  Dx even in [2, s-1],
                               so any value = 2 mod 4 in [2, 2s-4]
Write D0 = j (s-1) + b with b = D0 mod (s-1), which is EVEN since s-1 and
D0 are.  Then
   b = 0            -> j chords
   b = 2 mod 4      -> j + 2 chords   (one crossing pair)
   b = 0 mod 4, b>0 -> j + 4 chords   (two pairs, B = 2 and b-2)
so at most floor(D0/(s-1)) + 4 chords, which is O(1 + D0/s).
The only question is whether the pieces FIT on the circle.
"""
import sys
from dual import element, is_cycle
from skip4 import uniform

def bigpow(n):
    L = 1
    while L*2 <= n: L *= 2
    return L if L < n else n//2

def build(n, s, D0):
    """returns (chords, chord count) or None if it does not fit."""
    A = s-1
    j = D0 // A
    b = D0 - j*A
    pairs = []
    if b == 0:
        pass
    elif b % 4 == 2:
        pairs = [b]
    else:
        pairs = [2, b-2]
    if any(B < 2 or B > 2*s-4 for B in pairs): return None
    # lay out: j disjoint chords, then the pairs
    chords = []
    pos = 0
    for _ in range(j):
        chords.append((pos, pos+s)); pos += s+1
    for B in pairs:
        Dx = (B+2)//2
        if Dx % 2: return None
        chords.append((pos, pos+s))
        chords.append((pos+Dx, pos+Dx+s))
        pos += s+Dx+1
    if pos > n: return None
    out = []
    for (u, v) in chords:
        u %= n; v %= n
        out.append((min(u, v), max(u, v)))
    return out, len(chords)

if __name__ == '__main__':
    print("the construction, verified by building the cycle")
    print("   s    2N     L    D0    chords   built    length ok")
    tot = 0; okc = 0; nofit = []
    for s in range(5, 42, 2):
        for n in range(2*(s+2), 24*s+1, 2):
            L = bigpow(n); D0 = n-L
            tot += 1
            if D0 == 0: okc += 1; continue
            r = build(n, s, D0)
            if r is None: nofit.append((s, n, D0)); continue
            chords, k = r
            m = uniform(n, s)
            # every chord must really be a chord of U(n,s)
            if not all(m[u] == v for (u, v) in chords):
                nofit.append((s, n, D0, 'not a chord')); continue
            got = is_cycle(element(chords, n, 1), n)
            if got == L: okc += 1
            else: nofit.append((s, n, D0, 'len %s' % got))
    print("  pairs %d   construction gives length L on %d   failures %d"
          % (tot, okc, len(nofit)))
    print("  first failures:", nofit[:10])
