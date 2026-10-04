"""cl426.py -- closing (CL): the Triple Calculus, machine-verified.

For three even connecting paths P_0, P_1, P_2 (same terminal pair), write
the intersection-pattern parities:
    A = |P_0 ^ P_1 ^ P_2|, x = |P_0 ^ P_1 \\ P_2|, y = |P_0 ^ P_2 \\ P_1|,
    z = |P_1 ^ P_2 \\ P_0|,  p_i = private edges of P_i     (all mod 2).
With D_1 = P_0 xor P_1, D_2 = P_0 xor P_2:

  (T1)  lambda(D_1, D_2) = p_0 + z            (set identity, no SC needed)
  (T2)  path evenness:  p_0 = A + x + y,  p_1 = A + x + z,  p_2 = A + y + z
  (T3)  epsilon := lambda(D_1, D_2) = A + x + y + z
  (T4)  IF |P_0| == |P_1| == |P_2| (mod 4)  [the (SC) filter]:
            q(D_1) == 2(A + x),  q(D_2) == 2(A + y),
            q(D_1 + D_2) == 2(A + z)     (mod 4)
        so the failure plane's Brown type is (2,2,2) iff x == y == z, and
        the standard hyperbolic (2,0,0)-type iff exactly one agrees.
  (T5)  CLOSING triples (P_0 xor P_1 xor P_2 is again a path) have
        epsilon = 0 -- 222.2 recovered by the calculus.

Machine: sample even-path triples on random bridges; verify T1-T3 on all,
T4 on mod-4-aligned triples, T5 on closing triples; report the epsilon = 1
census and its Brown pattern split.
"""
import random

def machinery(VER, EDG):
    E = {frozenset(e): i for i, e in enumerate(EDG)}
    adj = {v: [] for v in VER}
    for u, w in EDG: adj[u].append(w); adj[w].append(u)
    bit = lambda u, w: 1 << E[frozenset((u, w))]
    return E, adj, bit

pop = lambda m: bin(m).count('1')

def all_paths(adj, bit, src, dst, cap=6000):
    out = []
    def dfs(v, mask, seen):
        if len(out) >= cap: return
        if v == dst: out.append(mask); return
        for w in adj[v]:
            if w not in seen: dfs(w, mask | bit(v, w), seen | {w})
    dfs(src, 0, {src})
    return out

def is_path(mask, E, src, dst):
    deg = {}
    for e, i in E.items():
        if mask >> i & 1:
            for v in e: deg[v] = deg.get(v, 0) + 1
    if not deg: return False
    odd = [v for v, d in deg.items() if d % 2 == 1]
    if sorted(map(str, odd)) != sorted(map(str, [src, dst])): return False
    if any(d != 2 for v, d in deg.items() if v not in (src, dst)): return False
    if deg.get(src) != 1 or deg.get(dst) != 1: return False
    # connectivity of the support
    verts = list(deg); vi = {v: i for i, v in enumerate(verts)}
    seen = {verts[0]}; st = [verts[0]]
    adj2 = {v: [] for v in verts}
    for e, i in E.items():
        if mask >> i & 1:
            a, b = tuple(e); adj2[a].append(b); adj2[b].append(a)
    while st:
        v = st.pop()
        for w in adj2[v]:
            if w not in seen: seen.add(w); st.append(w)
    return len(seen) == len(verts)

def random_bridge(n):
    for _ in range(4000):
        stubs = [v for v in range(n) for _ in range(3)]
        random.shuffle(stubs)
        E = set(); ok = True
        for i in range(0, len(stubs), 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b or frozenset((a, b)) in E: ok = False; break
            E.add(frozenset((a, b)))
        if not ok: continue
        adj = {v: [] for v in range(n)}
        for e in E:
            a, b = tuple(e); adj[a].append(b); adj[b].append(a)
        g = 99
        for s in range(n):
            dist = {s: 0}; par = {s: None}; q = [s]
            while q:
                v = q.pop(0)
                for w in adj[v]:
                    if w not in dist:
                        dist[w] = dist[v] + 1; par[w] = v; q.append(w)
                    elif par[v] != w: g = min(g, dist[v] + dist[w] + 1)
        if g < 5: continue
        def conn(rm):
            vs = [v for v in range(n) if v != rm]
            seen = {vs[0]}; st = [vs[0]]
            while st:
                v = st.pop()
                for w in adj[v]:
                    if w != rm and w not in seen: seen.add(w); st.append(w)
            return len(seen) == len(vs)
        if not all(conn(v) for v in range(n)): continue
        El = [tuple(e) for e in E]; random.shuffle(El)
        picked = []
        for e in El:
            if all(v not in f for f in picked for v in e): picked.append(e)
            if len(picked) == 3: break
        if len(picked) < 3: continue
        VER = list(range(n)) + ['t0', 't1', 't2']
        EDG = [tuple(e) for e in E
               if tuple(sorted(e)) not in [tuple(sorted(x)) for x in picked]]
        for k, (a, b) in enumerate(picked):
            EDG += [(a, 't%d' % k), ('t%d' % k, b)]
        return VER, EDG
    return None

random.seed(426)
stats = dict(triples=0, t123=0, sc=0, t4=0, close=0, t5=0,
             eps1=0, brown4=0, brown0=0)
graphs = 0
for n in [12, 14, 14, 16, 16, 18]:
    res = random_bridge(n)
    if not res: continue
    VER, EDG = res
    E, adj, bit = machinery(VER, EDG)
    graphs += 1
    for (s, t) in [('t0', 't1'), ('t0', 't2'), ('t1', 't2')]:
        ev = [m for m in all_paths(adj, bit, s, t) if pop(m) % 2 == 0]
        random.shuffle(ev)
        ev = ev[:40]
        for i in range(len(ev)):
            for j in range(i + 1, min(i + 8, len(ev))):
                for k in range(j + 1, min(j + 8, len(ev))):
                    P0, P1, P2 = ev[i], ev[j], ev[k]
                    stats['triples'] += 1
                    A = pop(P0 & P1 & P2)
                    x = pop(P0 & P1 & ~P2); y = pop(P0 & P2 & ~P1)
                    z = pop(P1 & P2 & ~P0)
                    p0 = pop(P0 & ~P1 & ~P2); p1 = pop(P1 & ~P0 & ~P2)
                    p2 = pop(P2 & ~P0 & ~P1)
                    D1, D2 = P0 ^ P1, P0 ^ P2
                    eps = pop(D1 & D2) % 2
                    ok123 = (eps == (p0 + z) % 2 == (A + x + y + z) % 2
                             and p0 % 2 == (A + x + y) % 2
                             and p1 % 2 == (A + x + z) % 2
                             and p2 % 2 == (A + y + z) % 2)
                    stats['t123'] += ok123
                    assert ok123
                    if pop(P0) % 4 == pop(P1) % 4 == pop(P2) % 4:
                        stats['sc'] += 1
                        q1, q2 = pop(D1) % 4, pop(D2) % 4
                        q3 = pop(D1 ^ D2) % 4
                        ok4 = (q1 == 2 * ((A + x) % 2)
                               and q2 == 2 * ((A + y) % 2)
                               and q3 == 2 * ((A + z) % 2))
                        stats['t4'] += ok4
                        assert ok4
                        if eps:
                            stats['eps1'] += 1
                            if (A + x) % 2 == (A + y) % 2 == (A + z) % 2:
                                stats['brown4'] += 1
                            else:
                                stats['brown0'] += 1
                        P3 = P0 ^ P1 ^ P2
                        if is_path(P3, E, s, t) and pop(P3) % 4 == pop(P0) % 4:
                            stats['close'] += 1
                            stats['t5'] += (eps == 0)
                            assert eps == 0
                        elif is_path(P3, E, s, t):
                            stats['close4'] = stats.get('close4', 0) + 1
print('bridges: %d, triples sampled: %d' % (graphs, stats['triples']))
print(' T1-T3 (set identities + evenness): %d/%d' %
      (stats['t123'], stats['triples']))
print(' mod-4-aligned triples: %d;  T4 (q-pattern formula): %d/%d' %
      (stats['sc'], stats['t4'], stats['sc']))
print(' fully aligned closing triples: %d;  T5 (=> eps = 0): %d/%d;'
      ' closing but 4th-path-misaligned (no claim): %d' %
      (stats['close'], stats['t5'], stats['close'], stats.get('close4', 0)))
print(' eps = 1 aligned planes: %d  (Brown-4: %d, Brown-0: %d)' %
      (stats['eps1'], stats['brown4'], stats['brown0']))
