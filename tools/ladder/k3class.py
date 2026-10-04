"""Classify the THREE-chord (M-2)-generators via the chain normal form.

A k=3 deficiency-2 element is: host minus three disjoint site arcs plus
three chords pairing the six site boundaries, no chord inside one site,
each site's two boundaries on different chords.  Gap profiles (3,1,1)
or (2,2,1) in some placement.  This enumerates, over site positions and
all boundary pairings, which (profile placement, pairing) combos give a
single (M-2)-cycle -- and checks the verdict is POSITION-FREE (depends
only on the combinatorial type, not where the sites sit).
"""
import sys
from itertools import permutations
from descend2 import cycles_of_len

def valid(M, starts, gaps, pairing):
    """sites at starts[i]..starts[i]+gaps[i]; pairing = perfect matching
    on boundary indices 0..5 (site i has boundaries 2i = left = starts[i],
    2i+1 = right = starts[i]+gaps[i]).  Returns True iff the element is a
    single simple (M-2)-cycle with girth-legal chords."""
    pos = []
    for i in range(3):
        pos.append(starts[i] % M)
        pos.append((starts[i] + gaps[i]) % M)
    # sites must be disjoint (with a vertex of slack between them so
    # boundaries are distinct and arcs do not touch)
    used = set()
    for i in range(3):
        for t in range(gaps[i] + 1):
            v = (starts[i] + t) % M
            if v in used:
                return None
            used.add(v)
    chords = []
    for a, b in pairing:
        u, v = pos[a], pos[b]
        if u == v:
            return None
        s = min((u - v) % M, (v - u) % M)
        if s < 4:
            return None                    # girth
        chords.append(tuple(sorted((u, v))))
    if len(set(chords)) < 3:
        return None
    # the element: host minus site edges plus chords -- test via the
    # generic detector, requiring all three chords used
    for e in cycles_of_len(M, sorted(chords), M - 2):
        es = set(tuple(sorted(t)) for t in e)
        if all(c in es for c in chords):
            # confirm the omitted vertices are the site interiors
            return True
    return False

def pairings():
    """perfect matchings on {0..5} with no pair inside one site and each
    site's boundaries on different chords (automatic given the first)."""
    out = []
    items = list(range(6))
    def rec(rem, cur):
        if not rem:
            out.append(tuple(cur)); return
        a = rem[0]
        for b in rem[1:]:
            if a // 2 == b // 2:
                continue                   # intra-site
            rec([x for x in rem[1:] if x != b], cur + [(a, b)])
    rec(items, [])
    return out

if __name__ == '__main__':
    M = 30
    PROFILES = sorted(set(permutations((3, 1, 1))) | set(permutations((2, 2, 1))))
    PAIRS = pairings()
    print('profiles:', PROFILES)
    print('pairings on 6 boundaries (inter-site):', len(PAIRS))
    verdicts = {}
    for prof in PROFILES:
        for pg in PAIRS:
            results = set()
            # sweep site positions: s1 = 0, s2, s3 varying with slack
            tried = 0
            for s2 in range(prof[0] + 2, M - prof[1] - prof[2] - 4, 1):
                for s3 in range(s2 + prof[1] + 2, M - prof[2] - 1, 1):
                    r = valid(M, (0, s2, s3), prof, pg)
                    if r is None:
                        continue
                    tried += 1
                    results.add(r)
                    if tried >= 40 and len(results) == 1:
                        break
                if tried >= 40 and len(results) == 1:
                    break
            if not results:
                continue
            verdicts[(prof, pg)] = results
    good = [(p, g) for (p, g), r in verdicts.items() if r == {True}]
    mixed = [(p, g) for (p, g), r in verdicts.items() if len(r) > 1]
    bad = [(p, g) for (p, g), r in verdicts.items() if r == {False}]
    print('valid generator types (always single cycle):', len(good))
    for p, g in sorted(good):
        print('   profile %s  pairing %s' % (p, g))
    print('never-valid types:', len(bad))
    print('POSITION-DEPENDENT types (mixed verdicts):', len(mixed))
    for p, g in sorted(mixed):
        print('   MIXED: profile %s  pairing %s' % (p, g))
