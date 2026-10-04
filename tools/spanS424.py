"""spanS424.py -- the proof machinery for (S), validated piecewise.

THEOREM (S): B 2-connected, non-bipartite, three degree-2 terminals ==>
W := span of even path-differences over the three pairs = C_0.

Proof architecture (424):
 L1  TWO-PATH LEMMA: every cycle Z is P1 xor P2 for a single pair (Menger
     connectors).  ==> V_ij := span of ALL-parity differences = C.
 L2  PAIR-COMPLETENESS: fixing even/odd bases (parity freedom, 326.2),
     V_ij = (W_ij + W'_ij) (+) <even xor odd base>, an odd-weight line, so
        W_ij + W'_ij = C_0   for EACH pair    (codimension argument).
 L3  FLIP CALCULUS: odd ij-paths carry signatures (arrival edge at b_i,
     arrival edge at b_j).  Moves absorbing W' into W:
       (i)   same-signature pairs: append a common ODD tail through the
             third terminal at either end -- odd.odd concatenation at a
             degree-2 terminal is an EVEN path of another pair;
       (ii)  single-end flips: tail at the AGREEING end;
       (iii) double crosses: the TRANSFER (R.A) xor (R'.A') moves the
             class to a single-end flip with free far-edge choice.
     All three reduce to the LINKING LEMMA (L): availability of the
     required disjoint odd tails/transfers, possibly after chaining.

Machine content here:
 (1) Pair-Completeness verified exactly, per pair, on random bridges.
 (2) Tail availability (move i/ii) and transfer availability (move iii)
     sampled on random bridges.
 (3) (S) itself re-verified per bridge.
"""
import random

def machinery(VER, EDG):
    E = {frozenset(e): i for i, e in enumerate(EDG)}
    adj = {v: [] for v in VER}
    for u, w in EDG: adj[u].append(w); adj[w].append(u)
    bit = lambda u, w: 1 << E[frozenset((u, w))]
    tree = set(); seen = {VER[0]}; st = [VER[0]]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w); tree.add(frozenset((v, w))); st.append(w)
    par = {VER[0]: None}; st = [VER[0]]
    while st:
        v = st.pop()
        for w in adj[v]:
            if frozenset((v, w)) in tree and w not in par:
                par[w] = v; st.append(w)
    def tp(a, b):
        m = 0; anc = set(); x = a
        while x is not None: anc.add(x); x = par[x]
        y = b
        while y not in anc: m |= bit(y, par[y]); y = par[y]
        x = a
        while x != y: m |= bit(x, par[x]); x = par[x]
        return m
    fund = [(1 << E[e]) | tp(*tuple(e)) for e in E if e not in tree]
    return E, adj, bit, fund

pop = lambda m: bin(m).count('1')
def rankv(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return basis
def inspan(v, basis):
    for b in sorted(basis, reverse=True): v = min(v, v ^ b)
    return v == 0

def all_paths(adj, bit, src, dst):
    out = []
    def dfs(v, mask, seen, first):
        if v == dst: out.append((mask, first, prev[0])); return
        for w in adj[v]:
            if w not in seen:
                prev0 = prev[0]; prev[0] = w if False else prev[0]
                dfs(w, mask | bit(v, w), seen | {w},
                    first if first else w)
    prev = [None]
    # simpler: record path vertex list
    out2 = []
    def dfs2(v, mask, seen, order):
        if v == dst: out2.append((mask, order[:])); return
        for w in adj[v]:
            if w not in seen:
                order.append(w)
                dfs2(w, mask | bit(v, w), seen | {w}, order)
                order.pop()
    dfs2(src, 0, {src}, [src] if False else [src])
    return out2

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
                    if w != rm and w not in seen:
                        seen.add(w); st.append(w)
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

random.seed(424)
PAIRS = [('t0', 't1', 't2'), ('t0', 't2', 't1'), ('t1', 't2', 't0')]
tot = {'pc': 0, 'pcok': 0, 'tail': 0, 'tailok': 0, 'tr': 0, 'trok': 0,
       'S': 0, 'Sok': 0}
graphs = 0
for n in [12, 14, 14, 16, 16, 18, 18, 16]:
    res = random_bridge(n)
    if not res: continue
    VER, EDG = res
    E, adj, bit, fund = machinery(VER, EDG)
    mu = len(fund)
    odd = [c for c in fund if pop(c) % 2 == 1]
    if not odd: continue
    graphs += 1
    C0 = rankv([c for c in fund if pop(c) % 2 == 0] +
               [a ^ b for a in odd for b in odd])
    P = {}
    for (s, t, k) in PAIRS:
        P[(s, t)] = all_paths(adj, bit, s, t)
    def getP(a, b):
        if (a, b) in P: return P[(a, b)]
        return [(m, o[::-1]) for m, o in P[(b, a)]]
    # (1) Pair-Completeness per pair
    for (s, t, k) in PAIRS:
        ps = P[(s, t)]
        ev = [m for m, o in ps if pop(m) % 2 == 0]
        od = [m for m, o in ps if pop(m) % 2 == 1]
        if not ev or not od: continue
        Wij = rankv([m ^ ev[0] for m in ev[1:]])
        Wp = rankv([m ^ od[0] for m in od[1:]])
        tot['pc'] += 1
        if len(rankv(Wij + Wp)) == len(C0): tot['pcok'] += 1
    # (2) tail availability for same-signature odd pairs + transfers
    for (s, t, k) in PAIRS:
        od = [(m, o) for m, o in P[(s, t)] if pop(m) % 2 == 1]
        odk = [(m, o) for m, o in getP(t, k) if pop(m) % 2 == 1]
        random.shuffle(od)
        for a in range(min(6, len(od))):
            for b in range(a + 1, min(6, len(od))):
                (m1, o1), (m2, o2) = od[a], od[b]
                if o1[-2] != o2[-2]: continue        # same edge at t
                tot['tail'] += 1
                V1 = set(o1) | set(o2)
                found = False
                for mk, ok_ in odk:
                    if ok_[1] == o1[-2]: continue    # must use other t-edge
                    if set(ok_[1:]) & (V1 - {t}): continue
                    found = True; break
                tot['tailok'] += found
        # transfers: double-cross pairs, far-edge-matched odd ki pair
        oki = [(m, o) for m, o in getP(k, s) if pop(m) % 2 == 1]
        cnt = 0
        for a in range(min(6, len(od))):
            for b in range(a + 1, min(6, len(od))):
                (m1, o1), (m2, o2) = od[a], od[b]
                if o1[1] == o2[1] or o1[-2] == o2[-2]: continue  # dbl cross
                if cnt >= 4: break
                cnt += 1; tot['tr'] += 1
                found = False
                for mR, oR in oki:
                    if oR[-2] == o1[1]: continue     # R arrives via other
                    if set(oR[:-1]) & set(o1): continue
                    for mR2, oR2 in oki:
                        if oR2[-2] == o2[1]: continue
                        if oR2[1] != oR[1]: continue  # same far edge at k
                        if set(oR2[:-1]) & set(o2): continue
                        found = True; break
                    if found: break
                tot['trok'] += found
    # (3) (S)
    diffs = []
    for (s, t, k) in PAIRS:
        ev = [m for m, o in P[(s, t)] if pop(m) % 2 == 0]
        if len(ev) > 1: diffs += [m ^ ev[0] for m in ev[1:]]
    tot['S'] += 1
    if len(rankv(diffs)) == len(C0): tot['Sok'] += 1

print('bridges: %d' % graphs)
print('L2 Pair-Completeness (W_ij + W\'_ij = C_0, per pair): %d/%d'
      % (tot['pcok'], tot['pc']))
print('move (i/ii) common odd tail available:  %d/%d sampled pairs'
      % (tot['tailok'], tot['tail']))
print('move (iii) transfer available:          %d/%d sampled crosses'
      % (tot['trok'], tot['tr']))
print('(S) W = C_0:                            %d/%d bridges'
      % (tot['Sok'], tot['S']))
