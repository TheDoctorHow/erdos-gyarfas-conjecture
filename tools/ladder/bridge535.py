import random
from block534 import (matchings, open_legal, girth, ladder_adj, spectrum,
                      crossing_components)
from path535 import side_graph, path_lengths, cycle_lengths_side

B12 = [m for m in matchings(list(range(1, 13))) if open_legal(m)]
B14 = [m for m in matchings(list(range(1, 15))) if open_legal(m)]
print("blocks: 12 ->", len(B12), "  14 ->", len(B14))
print("max span in a 12-block:", max(b - a for m in B12 for a, b in m))
print("max span in a 14-block:", max(b - a for m in B14 for a, b in m))
print("so cross constraint s_e+s_f <= 2N-3 at 2N=24 needs <=21, worst 9+9=18")
print("                              at 2N=26 needs <=23, worst 9+13=22")
print()

def assemble(m1, m2, twoN):
    k = 2 * len(m1)
    return sorted([(a, b) for (a, b) in m1] +
                  [(a + k, b + k) for (a, b) in m2])

def predicted(m1, m2, twoN):
    k = 2 * len(m1)
    g1 = side_graph(k, m1)
    g2 = side_graph(twoN - k, m2)
    P1 = path_lengths(g1, 1, k)
    P2 = path_lengths(g2, 1, twoN - k)
    S = set(cycle_lengths_side(g1, k)) | set(cycle_lengths_side(g2, twoN - k))
    for p in P1:
        for q in P2:
            S.add(2 + p + q)
    return S

print("BRIDGE DECOMPOSITION  spectrum = C_I u C_J u (2 + P_I + P_J)")
random.seed(11)
tests = [(B12[0], B12[0], 24), (B12[7], B12[19], 24),
         (B12[0], B14[0], 26), (B12[3], B14[500], 26)]
for _ in range(6):
    tests.append((random.choice(B12), random.choice(B12), 24))
    tests.append((random.choice(B12), random.choice(B14), 26))
ok = 0
for (m1, m2, twoN) in tests:
    ch = assemble(m1, m2, twoN)
    adj = ladder_adj(twoN, ch)
    brute = set(spectrum(adj))
    pred = predicted(m1, m2, twoN)
    same = (brute == pred)
    ok += same
    if not same:
        print("   MISMATCH", twoN, sorted(brute ^ pred))
print("   agree on %d of %d assemblies (brute force vs formula)" %
      (ok, len(tests)))
print()

# ---- (LQ3) over the WHOLE disconnected class at 2N = 24 and 26 ----
def descents(m2, ch):
    return set(m2 - 1 - p for p in path_lengths(side_graph(m2, ch), 1, m2))

D12 = [descents(12, m) for m in B12]
D14 = [descents(14, m) for m in B14]
C12 = [set(cycle_lengths_side(side_graph(12, m), 12)) for m in B12]

def has_pow2(twoN, Di, Dj, Ci, Cj):
    for t in range(3, 12):
        p = 1 << t
        if p > twoN:
            break
        if p in Ci or p in Cj:
            return True
        need = twoN - p
        if any(need - x in Dj for x in Di):
            return True
    return False

f24 = sum(1 for i in range(len(B12)) for j in range(len(B12))
          if not has_pow2(24, D12[i], D12[j], C12[i], C12[j]))
print("2N = 24: shape pairs failing (LQ3):", f24, "of", len(B12) ** 2)
print("   whole class = 12 placements x 34 x 34 =", 12 * 34 * 34, "ladders")

C14 = [set(cycle_lengths_side(side_graph(14, m), 14)) for m in B14]
f26 = sum(1 for i in range(len(B12)) for j in range(len(B14))
          if not has_pow2(26, D12[i], D14[j], C12[i], C14[j]))
print("2N = 26: shape pairs failing (LQ3):", f26, "of",
      len(B12) * len(B14))
print("   whole class = 26 placements x 34 x 792 =", 26 * 34 * 792,
      "ladders")
print()

# independent confirmation on a random sample, by brute-force spectrum
random.seed(7)
bad = 0
for _ in range(40):
    twoN = random.choice([24, 26])
    m1 = random.choice(B12)
    m2 = random.choice(B12 if twoN == 24 else B14)
    ch = assemble(m1, m2, twoN)
    adj = ladder_adj(twoN, ch)
    if girth(adj) < 5:
        print("   ILLEGAL assembly!", twoN, ch)
        bad += 1
        continue
    sp = set(spectrum(adj))
    if not any(x & (x - 1) == 0 for x in sp):
        print("   NO POWER OF TWO", twoN, ch)
        bad += 1
print("random sample of 40 assemblies, brute-force girth+spectrum:",
      "all legal and all with a power of two" if bad == 0 else "%d bad" % bad)

# which powers of two, and how often the block alone suffices
alone = sum(1 for i in range(len(B12)) if 8 in C12[i])
print("12-blocks whose own side holds an 8-cycle: %d of %d" % (alone, len(B12)))
alone14 = sum(1 for j in range(len(B14))
              if any(p in C14[j] for p in (8, 16)))
print("14-blocks whose own side holds an 8- or 16-cycle: %d of %d"
      % (alone14, len(B14)))
