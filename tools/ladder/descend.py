"""Test the DEFICIENCY-2 TELESCOPE.

A cycle of length exactly L on 2N vertices omits exactly 2N - L
vertices.  So the sharpened statement telescopes: find a cycle omitting
TWO vertices, re-present the graph with that cycle as the new host (575:
the presentation is a choice, not a property), and repeat; at size L the
host circle itself is the target.  Girth is inherited at every step
because each new object is a SUBGRAPH of the old one - which is exactly
what 590's contraction induction lacked.  The costs: each step loses the
chords dangling into the two omitted vertices, so the partial matching
thins as the telescope descends.

Measured here on random girth-legal ladders:
  A. single step: how many (M-2)-cycles there are, and what fraction of
     the induced partial ladders H retain an L-cycle;
  B. full descent: whether SOME chain of deficiency-2 steps reaches the
     band bottom (DFS, first success);
  C. the deficit trajectory: chords at each level on the successful path.
"""
import sys, time, random
from pushL import bigpow
from hitset import random_ladder

def elements(M, ch):
    """every simple cycle of host-circle-M + chords, as (length, edges);
    (S,eps) scan over the full cycle space, 2^(k+1) elements."""
    k = len(ch)
    out = []
    for mask in range(1 << k):
        S = [ch[i] for i in range(k) if mask >> i & 1]
        tog = [0] * M
        for a, b in S:
            tog[a] ^= 1; tog[b] ^= 1
        for eps in (0, 1):
            t = eps
            edges = list(S)
            for v in range(M):
                t ^= tog[v]
                if t:
                    edges.append((v, (v + 1) % M))
            if not edges:
                continue
            deg = {}
            for a, b in edges:
                deg[a] = deg.get(a, 0) + 1
                deg[b] = deg.get(b, 0) + 1
            if any(d != 2 for d in deg.values()):
                continue
            adj = {}
            for a, b in edges:
                adj.setdefault(a, []).append(b)
                adj.setdefault(b, []).append(a)
            st = [next(iter(adj))]; seen = {st[0]}
            while st:
                u = st.pop()
                for w in adj[u]:
                    if w not in seen:
                        seen.add(w); st.append(w)
            if len(seen) != len(deg):
                continue
            out.append((len(edges), edges))
    return out

def cycles_of_len(M, ch, target):
    return [e for ln, e in elements(M, ch) if ln == target]

def has_len(M, ch, target):
    return any(ln == target for ln, e in elements(M, ch))

def step(M, ch, cyc_edges):
    """re-present: cyc_edges (length M-2) becomes the new host; return
    (M-2, new chords).  New object is a subgraph, so girth is kept."""
    es = set(tuple(sorted(e)) for e in cyc_edges)
    verts = set()
    for a, b in es:
        verts.add(a); verts.add(b)
    adjfull = {v: set() for v in range(M)}
    for v in range(M):
        adjfull[v].add((v + 1) % M); adjfull[(v + 1) % M].add(v)
    for a, b in ch:
        adjfull[a].add(b); adjfull[b].add(a)
    cyc_adj = {}
    for a, b in es:
        cyc_adj.setdefault(a, []).append(b)
        cyc_adj.setdefault(b, []).append(a)
    start = min(verts)
    order = [start]
    prev, cur = None, start
    while True:
        nxt = next(w for w in cyc_adj[cur] if w != prev)
        if nxt == start:
            break
        order.append(nxt); prev, cur = cur, nxt
    assert len(order) == len(verts)
    pos = {v: i for i, v in enumerate(order)}
    newch = set()
    for v in order:
        for w in adjfull[v]:
            e = tuple(sorted((v, w)))
            if e in es:
                continue
            if w in pos:
                newch.add(tuple(sorted((pos[v], pos[w]))))
    return len(order), sorted(newch)

def telescope(M, ch, L, trail, node_budget):
    """DFS: does some deficiency-2 chain reach the band bottom?"""
    if M == L:
        return True
    if node_budget[0] <= 0:
        return False
    node_budget[0] -= 1
    cands = cycles_of_len(M, ch, M - 2)
    trail.append((M, len(ch), len(cands)))
    for e in cands:
        Mn, nch = step(M, ch, e)
        if telescope(Mn, nch, L, trail, node_budget):
            return True
    return False

if __name__ == '__main__':
    mode = sys.argv[1]
    rng = random.Random(int(sys.argv[2]))
    trials = int(sys.argv[3])
    n = int(sys.argv[4])
    L = bigpow(n)
    if mode == 'A':
        # single step at size n: (n-2)-cycles, and retention of L in H
        tot = good = 0
        worst = None
        for t in range(trials):
            ch = random_ladder(n, rng)
            if ch is None:
                continue
            cands = cycles_of_len(n, ch, n - 2)
            keep = 0
            for e in cands:
                Mn, nch = step(n, ch, e)
                if has_len(Mn, nch, L):
                    keep += 1
            tot += 1
            if cands and keep:
                good += 1
            frac = (keep / len(cands)) if cands else 0.0
            if worst is None or frac < worst[0]:
                worst = (frac, len(cands), keep)
            print("  ladder %2d: %3d (n-2)-cycles, %3d descend well (%.2f)"
                  % (t, len(cands), keep, frac))
            sys.stdout.flush()
        print("  == %d/%d ladders have a good step; worst fraction %.2f "
              "(%d of %d)" % (good, tot, worst[0], worst[2], worst[1]))
    elif mode == 'B':
        ok = 0; tot = 0
        for t in range(trials):
            ch = random_ladder(n, rng)
            if ch is None:
                continue
            trail = []
            budget = [4000]
            r = telescope(n, ch, L, trail, budget)
            tot += 1; ok += bool(r)
            path = trail[:6]
            print("  ladder %2d: %s  nodes=%d  first-levels %s"
                  % (t, "REACHED %d" % L if r else "FAILED",
                     4000 - budget[0], path))
            sys.stdout.flush()
        print("  == telescope succeeded on %d/%d ladders" % (ok, tot))
