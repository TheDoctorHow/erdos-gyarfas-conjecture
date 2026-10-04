"""bip294.py -- test the SHARPENED form of (E-CORE).

  (BIP)  a terminal pair (b_i,b_j) FAILS (L2**)
         <=>  some edge e has B - e BIPARTITE, with b_i and b_j in DIFFERENT classes.

The <= direction is a two-line proof (see below).  This tests =>, and the consequence
that matters: since three terminals cannot be pairwise bichromatic in a 2-colouring,
AT MOST TWO of the three pairs can fail -- so some pair satisfies (L2**).

Also measured: tau_odd^edge(B) = min edges whose deletion makes B bipartite, and the
number of failing pairs per bridge (predicted 0 or 2, never 1 or 3).
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split('"""', 2)[2])

def bipartition(a):
    """returns colouring dict or None"""
    col = {}
    for s in a:
        if s in col: continue
        col[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in col: col[w] = 1 - col[u]; q.append(w)
                elif col[w] == col[u]: return None
    return col

def tau_edge(E, cap=3):
    """min number of edges to delete to make the graph bipartite (up to cap)"""
    El = sorted(E, key=lambda s: sorted(map(str, s)))
    for k in range(cap + 1):
        for S in itertools.combinations(El, k):
            if bipartition(adj_of(set(El) - set(S))) is not None:
                return k, set(S)
    return None, None

def even_paths_all(a, s, t):
    out = []; path = [s]; onp = {s}
    def rec(u):
        for w in a[u]:
            if w == t:
                p = path + [t]
                if (len(p) - 1) % 2 == 0: out.append(p)
            elif w not in onp:
                onp.add(w); path.append(w); rec(w); path.pop(); onp.discard(w)
    rec(s); return out

def main():
    seed = int(sys.argv[1]); mlist = [int(x) for x in sys.argv[2].split(',')]
    per = int(sys.argv[3])
    rng = random.Random(seed)
    st = collections.Counter()
    perbridge = collections.Counter()
    for m in mlist:
        got = 0; att = 0
        while got < per and att < per * 200:
            att += 1
            a, terms = make_bridge(m, rng)
            if a is None: continue
            got += 1; st['bridges'] += 1
            E = {frozenset((u, w)) for u in a for w in a[u]}
            k, S = tau_edge(E)
            st['tau_edge_%s' % k] += 1
            # the bipartizing edges: all e with B-e bipartite
            bipedges = [e for e in E if bipartition(adj_of(E - {e})) is not None]
            st['num_bipartizing_edges_%d' % len(bipedges)] += 1
            nfail = 0
            for b1, b2 in itertools.combinations(terms, 2):
                ev = even_paths_all(a, b1, b2)
                if len(ev) < 2: st['degenerate'] += 1; continue
                Es = [frozenset(frozenset((p[i], p[i+1])) for i in range(len(p)-1))
                      for p in ev]
                core = set(Es[0])
                for x in Es[1:]: core &= x
                if core:
                    fail = True
                else:
                    fail = not any(not (Es[i] & Es[j])
                                   for i in range(len(Es)) for j in range(i+1, len(Es)))
                st['pairs'] += 1
                if not fail:
                    st['HOLDS'] += 1
                    # (BIP) <= direction must never be violated
                    for e in bipedges:
                        c = bipartition(adj_of(E - {e}))
                        if c[b1] != c[b2]:
                            st['*** BIP EASY DIRECTION VIOLATED ***'] += 1
                    continue
                nfail += 1; st['FAIL'] += 1
                witness = [e for e in bipedges
                           if bipartition(adj_of(E - {e}))[b1]
                           != bipartition(adj_of(E - {e}))[b2]]
                if witness:
                    st['fail_with_BIP_witness'] += 1
                    st['witness_is_the_core_edge'] += (set(witness) == core)
                else:
                    st['*** BIP COUNTEREXAMPLE ***'] += 1
                    print('BIP COUNTEREXAMPLE n=%d pair=%s,%s core=%s bipedges=%s'
                          % (len(a), b1, b2,
                             sorted(tuple(sorted(map(str, x))) for x in core),
                             sorted(tuple(sorted(map(str, x))) for x in bipedges)))
                    print('   graph:', sorted(tuple(sorted(map(str, x))) for x in E))
            perbridge[nfail] += 1
    print('m=%s seed=%d' % (mlist, seed))
    for k in sorted(st): print('  %-34s %s' % (k, st[k]))
    print('  failing pairs per bridge:', dict(sorted(perbridge.items())))

main()
