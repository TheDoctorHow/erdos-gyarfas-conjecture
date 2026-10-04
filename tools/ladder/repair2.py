"""Stage 2: seed from near-miss designs, perturb locally, larger search."""
import random, subprocess, sys
from repair import (iv, N, fixed, pools, mixed, odd_feet, even_feet,
                    corridors_ok, to_g6, rand_matching, build)

rng = random.Random(99)

def matching_ok(pairs):
    return all(corridors_ok(a, b) for a, b in pairs)

def perturb(pairs, rng, k=2):
    """Re-pair k randomly chosen pairs among themselves."""
    for _ in range(100):
        p = [list(t) for t in pairs]
        idx = rng.sample(range(len(p)), k)
        feet = [f for i in idx for f in p[i]]
        rng.shuffle(feet)
        newp = []
        ok = True
        f = feet[:]
        while f:
            a = f.pop()
            c = [i for i, b in enumerate(f) if corridors_ok(tuple(a), tuple(b))]
            if not c:
                ok = False
                break
            newp.append((tuple(a), tuple(f.pop(rng.choice(c)))))
        if not ok:
            continue
        out = [tuple(map(tuple, p[i])) for i in range(len(p)) if i not in idx]
        out += newp
        return out
    return None

# stage 1 rebuild: regenerate the same batch to recover near-miss designs
rng1 = random.Random(20260803)
seeds = []
lines = []
metas = []
for t in range(400):
    op = rand_matching(odd_feet, rng1)
    ep = rand_matching(even_feet, rng1)
    if op is None or ep is None:
        continue
    lines.append(to_g6(build(op, ep)))
    metas.append((op, ep))
with open("b1.g6", "w") as f:
    f.write("\n".join(lines) + "\n")
r = subprocess.run(["./egtool.exe", "-q"], stdin=open("b1.g6"),
                   capture_output=True, text=True)
near = set()
for ln in r.stdout.splitlines():
    if ln.startswith("NEARMISS"):
        near.add(ln.split()[1])
for g6, meta in zip(lines, metas):
    if g6 in near:
        seeds.append(meta)
print("near-miss seeds:", len(seeds))

# stage 2: perturb seeds
best = None
lines2, metas2 = [], []
for (op, ep) in seeds:
    for _ in range(120):
        which = rng.random()
        op2 = perturb(op, rng, rng.choice([2, 3])) if which < 0.5 else op
        ep2 = perturb(ep, rng, rng.choice([2, 3])) if which >= 0.5 else ep
        if op2 is None or ep2 is None:
            continue
        lines2.append(to_g6(build(list(op2), list(ep2))))
        metas2.append((op2, ep2))
print("stage-2 candidates:", len(lines2))
with open("b2.g6", "w") as f:
    f.write("\n".join(lines2) + "\n")
r = subprocess.run(["./egtool.exe", "-q", "-b", "2000000000"],
                   stdin=open("b2.g6"), capture_output=True, text=True)
print(r.stderr.strip())
hits = [ln for ln in r.stdout.splitlines() if ln.startswith("COUNTEREXAMPLE")]
nears2 = [ln for ln in r.stdout.splitlines() if ln.startswith("NEARMISS")]
print("hits:", len(hits), " nearmisses:", len(nears2))
if hits:
    g6hit = hits[0].split()[1]
    idx = lines2.index(g6hit)
    with open("hit96.g6", "w") as f:
        f.write(g6hit + "\n")
    print(hits[0])
    print("ODD :", metas2[idx][0])
    print("EVEN:", metas2[idx][1])
