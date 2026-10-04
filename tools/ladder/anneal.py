"""Greedy 2-opt annealing on the attachment matching.

State: perfect matching on the 90 interiors (45 attachment edges) over the
fixed skeleton. Score: number of cycles of length 4, 8, 16, 32, 64
(counted exactly via egcore-style DFS in C? no -- via egtool -c parse).
Moves: pick two attachment edges, re-pair their four feet the other legal way.
Greedy descent with random restarts from near-miss seeds.
"""
import random, subprocess, sys

N = 96

def iv(x, y, k):
    return 6 + ((x - 1) * 3 + (y - 1)) * 10 + (k - 1)

skeleton = []
for x in (1, 2, 3):
    for y in (1, 2, 3):
        skeleton.append((x - 1, iv(x, y, 1)))
        for k in range(1, 10):
            skeleton.append((iv(x, y, k), iv(x, y, k + 1)))
        skeleton.append((iv(x, y, 10), 3 + y - 1))

def foot(v):
    """global id -> (x, y, u-offset) or None for branch vertices"""
    if v < 6:
        return None
    v -= 6
    c, k = divmod(v, 10)
    x, y = divmod(c, 3)
    return (x + 1, y + 1, k + 1)

def to_g6(edges):
    adj = [[0] * N for _ in range(N)]
    for a, b in edges:
        adj[a][b] = adj[b][a] = 1
    bits = []
    for j in range(1, N):
        for i in range(j):
            bits.append(adj[i][j])
    while len(bits) % 6:
        bits.append(0)
    out = [chr(126), chr(((N >> 12) & 63) + 63), chr(((N >> 6) & 63) + 63),
           chr((N & 63) + 63)]
    for k2 in range(0, len(bits), 6):
        v = 0
        for b in bits[k2:k2 + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return "".join(out)

def score_batch(matchings):
    """Return list of (#bad-length cycle classes present, mask) per matching.
    Uses egtool -q output; score = number of 2^k lengths present (we want 0);
    tie-break not available cheaply, so score = count of k with cycle."""
    lines = [to_g6(skeleton + m) for m in matchings]
    fn = "an.g6"
    with open(fn, "w") as f:
        f.write("\n".join(lines) + "\n")
    r = subprocess.run(["./egtool.exe"], stdin=open(fn),
                       capture_output=True, text=True)
    scores = []
    for ln in r.stdout.splitlines():
        if " pow2=" not in ln:
            continue
        part = ln.split(" pow2=")[1].split()[0]
        ks = [p for p in part.split(",") if p]
        scores.append(len(ks))
    return scores

def legal_pair(a, b):
    fa, fb = foot(a), foot(b)
    (x1, y1, k1), (x2, y2, k2) = fa, fb
    if x1 == x2 and y1 == y2:
        d = abs(k1 - k2)
        return d >= 5 and d != 7          # d=7 -> C8 (1+d)
    if x1 != x2 and y1 != y2:
        # disjoint: parity law kills all 1-edge C16/C32/C64
        # (odd |a-b| in {7,9} and odd a+b in {9,13,19} hits)
        return (k1 - k2) % 2 == 0
    if x1 == x2:
        i, j = k1, k2
    else:
        i, j = 11 - k1, 11 - k2
    s = i + j
    if s in (3, 7, 15):                   # through-cycle 4, 8, 16
        return False
    return abs(i - j) != 2                # far-route 34+(i-j) -> C32

def move(m, rng):
    """2-opt: re-pair two random edges."""
    m2 = m[:]
    i, j = rng.sample(range(len(m2)), 2)
    (a, b), (c, d) = m2[i], m2[j]
    opts = []
    if legal_pair(a, c) and legal_pair(b, d):
        opts.append(((a, c), (b, d)))
    if legal_pair(a, d) and legal_pair(b, c):
        opts.append(((a, d), (b, c)))
    if not opts:
        return None
    e1, e2 = rng.choice(opts)
    m2[i], m2[j] = e1, e2
    return m2

def from_g6(s):
    p = 0
    assert ord(s[0]) == 126
    n = 0
    for i in range(1, 4):
        n = (n << 6) | (ord(s[i]) - 63)
    p = 4
    bitpos = 0
    edges = []
    for j in range(1, n):
        for i in range(j):
            byte = bitpos // 6
            bit = 5 - (bitpos % 6)
            if ((ord(s[p + byte]) - 63) >> bit) & 1:
                edges.append((i, j))
            bitpos += 1
    return edges

_seeds = None
def seed_pool():
    """Near-miss designs from the stage-2 batch (b2.g6)."""
    global _seeds
    if _seeds is not None:
        return _seeds
    r = subprocess.run(["./egtool.exe", "-q"], stdin=open("b2.g6"),
                       capture_output=True, text=True)
    near = []
    for ln in r.stdout.splitlines():
        if ln.startswith("NEARMISS"):
            near.append(ln.split()[1])
    _seeds = near
    print("seed pool: %d near-miss designs" % len(near))
    return _seeds

def seed_from_tsearch(rng):
    pool = seed_pool()
    if not pool:
        return None
    g6 = rng.choice(pool)
    sk = set((min(a, b), max(a, b)) for a, b in skeleton)
    att = [e for e in from_g6(g6)
           if (min(e[0], e[1]), max(e[0], e[1])) not in sk]
    assert len(att) == 45, len(att)
    return att

import os

def seed_from_ckpt(rng):
    """Resume from the best checkpoint on disk, if any."""
    cks = sorted(f for f in os.listdir(".") if f.startswith("ckpt-"))
    if not cks:
        return None
    with open(cks[0]) as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    g6 = rng.choice(lines)
    sk = set((min(a, b), max(a, b)) for a, b in skeleton)
    att = [e for e in from_g6(g6)
           if (min(e[0], e[1]), max(e[0], e[1])) not in sk]
    return att if len(att) == 45 else None

def checkpoint(m, s):
    """Append design to ckpt-<score>.g6 the moment it is reached."""
    with open("ckpt-%d.g6" % s, "a") as f:
        f.write(to_g6(skeleton + m) + "\n")

def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
    best_overall = None
    for restart in range(20):
        m = None
        if restart % 3 == 1:
            m = seed_from_ckpt(rng)   # every 3rd restart resumes best-so-far
        if m is None:
            m = seed_from_tsearch(rng)
        if m is None:
            continue
        s = score_batch([m])[0]
        print("restart %d: seed score %d" % (restart, s))
        sys.stdout.flush()
        stall = 0
        while s > 0 and stall < 80:
            cands = []
            for _ in range(32):
                m2 = move(m, rng)
                if m2 is not None:
                    cands.append(m2)
            if not cands:
                break
            ss = score_batch(cands)
            bi = min(range(len(ss)), key=lambda t: ss[t])
            if ss[bi] <= s:
                if ss[bi] < s:
                    stall = 0
                    print("  improved: %d -> %d" % (s, ss[bi]))
                    sys.stdout.flush()
                    checkpoint(cands[bi], ss[bi])
                else:
                    stall += 1
                m, s = cands[bi], ss[bi]
            else:
                stall += 1
        print("restart %d final score: %d" % (restart, s))
        sys.stdout.flush()
        if s == 0:
            g6 = to_g6(skeleton + m)
            with open("HIT.g6", "w") as f:
                f.write(g6 + "\n")
            print("ZERO SCORE -- candidate written to HIT.g6")
            print(g6)
            return
        if best_overall is None or s < best_overall[0]:
            best_overall = (s, m)
            checkpoint(m, s)
    if best_overall:
        print("best score:", best_overall[0])

if __name__ == "__main__":
    main()
