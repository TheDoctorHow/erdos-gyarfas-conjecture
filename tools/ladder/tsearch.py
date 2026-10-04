"""Joint template + long-range pairing search for the m=1 K3,3 candidate.

Per branch vertex: a template = 5 cross-corridor hatch edges on items
{A,B,C}x{1..5} (local frame) + 5 exports, satisfying local 1-edge and
2-edge legality. Six vertices get independent random templates.
Mixed edge placed among exports (pins: shared-frame {0,1} mod 4, sum=1 mod 8).
Remaining exports paired long-range. Oracle: egtool.exe.
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

# pool corridor rosters (A, B, C) per vertex; local->global offset maps
POOLS = [
    ("u", [(1,1), (1,2), (1,3)]),
    ("u", [(2,2), (2,3), (2,1)]),
    ("u", [(3,3), (3,1), (3,2)]),
    ("v", [(1,1), (2,1), (3,1)]),
    ("v", [(2,2), (3,2), (1,2)]),
    ("v", [(3,3), (1,3), (2,3)]),
]
def glob_off(side, k):
    return k if side == "u" else 11 - k

def tmpl_legal(edges):
    """edges: list of ((ci, i), (cj, j)) with ci != cj, local offsets."""
    for (a, i), (b, j) in edges:
        s = i + j
        if s in (3, 7):
            return False
    for e1 in range(len(edges)):
        for e2 in range(e1 + 1, len(edges)):
            (a1, i1), (b1, j1) = edges[e1]
            (a2, i2), (b2, j2) = edges[e2]
            f1 = {(a1): i1, (b1): j1}
            f2 = {(a2): i2, (b2): j2}
            shared = set(f1) & set(f2)
            if len(shared) == 2:
                d = abs(f1[a1] - f2[a1]) + abs(f1[b1] - f2[b1])
                if d in (2, 6):
                    return False
            elif len(shared) == 1:
                c = shared.pop()
                d = abs(f1[c] - f2[c])
                rest = [v for k2, v in f1.items() if k2 != c] + \
                       [v for k2, v in f2.items() if k2 != c]
                if d + sum(rest) in (2, 6):
                    return False
    return True

def rand_template(rng):
    items = [(c, k) for c in range(3) for k in range(1, 6)]
    for _ in range(300):
        it = items[:]
        rng.shuffle(it)
        edges = []
        ok = True
        pool = it[:]
        # pick 5 edges greedily, leave 5 exports
        for _ in range(5):
            placed = False
            rng.shuffle(pool)
            for ai in range(len(pool)):
                for bi in range(ai + 1, len(pool)):
                    a, b = pool[ai], pool[bi]
                    if a[0] == b[0]:
                        continue
                    cand = edges + [(a, b)]
                    if tmpl_legal(cand):
                        edges = cand
                        pool = [p for pi, p in enumerate(pool)
                                if pi not in (ai, bi)]
                        placed = True
                        break
                if placed:
                    break
            if not placed:
                ok = False
                break
        if ok:
            return edges, pool  # pool = 5 exports (corridor-idx, local offset)
    return None, None

def corridors_ok(f1, f2):
    (x1, y1, k1), (x2, y2, k2) = f1, f2
    if x1 == x2 and y1 == y2:
        return abs(k1 - k2) in (6, 8, 10)
    if x1 != x2 and y1 != y2:
        return True
    if x1 == x2:
        i, j = k1, k2
    else:
        i, j = 11 - k1, 11 - k2
    s = i + j
    return s != 3 and s % 8 != 7

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
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return "".join(out)

def gen_design(rng):
    edges = list(skeleton)
    exports = []  # (x, y, u-offset)
    for (side, ros) in POOLS:
        te, ex = rand_template(rng)
        if te is None:
            return None
        for (ci, i), (cj, j) in te:
            (x1, y1), (x2, y2) = ros[ci], ros[cj]
            edges.append((iv(x1, y1, glob_off(side, i)),
                          iv(x2, y2, glob_off(side, j))))
        for (ci, k) in ex:
            x, y = ros[ci]
            exports.append((x, y, glob_off(side, k)))
    # mixed edge: shared-vertex pair of exports, pins {0,1} mod 4, sum=1 mod 8
    odd = [f for f in exports if f[2] % 2 == 1]
    even = [f for f in exports if f[2] % 2 == 0]
    mixed = None
    cand = []
    for a in odd:
        for b in even:
            if a[0] == b[0] and a[1] != b[1]:      # share u
                i, j = a[2], b[2]
            elif a[1] == b[1] and a[0] != b[0]:    # share v
                i, j = 11 - a[2], 11 - b[2]
            else:
                continue
            if {i % 4, j % 4} == {0, 1} and (i + j) % 8 == 1 and i + j != 3:
                cand.append((a, b))
    if not cand:
        return None
    mixed = rng.choice(cand)
    odd.remove(mixed[0])
    even.remove(mixed[1])
    edges.append((iv(*mixed[0]), iv(*mixed[1])))
    # long-range pairing
    for feet in (odd, even):
        f = feet[:]
        for _ in range(100):
            rng.shuffle(f)
            pairs, rest, ok = [], f[:], True
            while rest:
                a = rest.pop()
                c = [i for i, b in enumerate(rest) if corridors_ok(a, b)]
                if not c:
                    ok = False
                    break
                pairs.append((a, rest.pop(rng.choice(c))))
            if ok:
                for a, b in pairs:
                    edges.append((iv(*a), iv(*b)))
                break
        else:
            return None
    # simplicity + cubic check
    es = set()
    deg = [0] * N
    for a, b in edges:
        e = (min(a, b), max(a, b))
        if e in es:
            return None
        es.add(e)
        deg[a] += 1
        deg[b] += 1
    if any(d != 3 for d in deg):
        return None
    return sorted(es)

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
ROUNDS = 6
BATCH = 500
for rd in range(ROUNDS):
    lines = []
    while len(lines) < BATCH:
        d = gen_design(rng)
        if d is not None:
            lines.append(to_g6(d))
    fn = "t%d.g6" % rd
    with open(fn, "w") as f:
        f.write("\n".join(lines) + "\n")
    r = subprocess.run(["./egtool.exe", "-q"], stdin=open(fn),
                       capture_output=True, text=True)
    hits = [ln for ln in r.stdout.splitlines() if ln.startswith("COUNTEREXAMPLE")]
    nears = [ln for ln in r.stdout.splitlines() if ln.startswith("NEARMISS")]
    print("round %d: %s | hits=%d nearmiss=%d" %
          (rd, r.stderr.splitlines()[0], len(hits), len(nears)))
    sys.stdout.flush()
    if hits:
        with open("HIT.g6", "w") as f:
            f.write(hits[0].split()[1] + "\n")
        print(hits[0])
        break
