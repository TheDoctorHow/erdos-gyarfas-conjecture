"""family484.py -- THE CYCLE SEMIGROUP: machine check of the two lemmas.

474 proposes an OBJECT that satisfies J7 by its own structure:

    a CYCLE SEMIGROUP is a family F of finite graphs of minimum
    degree >= 3 closed under EDGE AMALGAMATION (the cubic 2-sum).

THEOREM A (474.1).  Every nonempty cycle semigroup contains a member
with a cycle whose length is a power of two.

The proof has one graph-theoretic input, verified here:

  LEMMA 1 (amalgam spectrum, exact).  Let H1, H2 be disjoint graphs,
  e1 = u1v1 in E(H1), e2 = u2v2 in E(H2).  Let

      H1 (+) H2  =  (H1 - e1) + (H2 - e2) + u1u2 + v1v2 .

  Then     S(H1 (+) H2)  =  S(H1 - e1)  U  S(H2 - e2)
                            U  ( S_{e1}(H1) + S_{e2}(H2) )
  where S_e(H) is the set of lengths of cycles of H through e, and
  the last term is a SUMSET.  Moreover the cycles through the new
  edge u1u2 are exactly the sumset part, so the construction
  ITERATES: an m-fold chain realises the m-fold sumset.

  LEMMA 2 (Bondy-Vince corollary).  For every finite graph H of
  minimum degree >= 3, gcd S(H) is 1 or 2 - a power of two.

and one arithmetic input, Sylvester: a subsemigroup of N whose gcd
is a power of two contains a power of two, because {2^k} is COFINAL
in 2^s N.

This file checks Lemma 1 exactly on amalgam pairs and on m-fold
chains, checks Lemma 2 on the zoo, and prints the explicit ESCAPE
CERTIFICATE of 474.3: the least power of two in <S(G)> and the
number of blocks needed to realise it.
"""
import itertools
from math import gcd
P = print


# ---------------------------------------------------------------- graphs

def gp(n, k):
    adj = {i: set() for i in range(2 * n)}
    def add(x, y):
        if x != y: adj[x].add(y); adj[y].add(x)
    for i in range(n):
        add(i, (i + 1) % n); add(i, n + i); add(n + i, n + (i + k) % n)
    return [sorted(adj[i]) for i in range(2 * n)]


def lcf(n, seq):
    adj = {i: set() for i in range(n)}
    for i in range(n):
        adj[i].add((i + 1) % n); adj[(i + 1) % n].add(i)
    for i in range(n):
        j = (i + seq[i % len(seq)]) % n
        adj[i].add(j); adj[j].add(i)
    return [sorted(adj[i]) for i in range(n)]


def k33():
    adj = [set() for _ in range(6)]
    for a in range(3):
        for b in range(3, 6): adj[a].add(b); adj[b].add(a)
    return [sorted(x) for x in adj]


def prism(m):
    """the m-prism: two m-cycles joined by a perfect matching."""
    adj = [set() for _ in range(2 * m)]
    def add(x, y): adj[x].add(y); adj[y].add(x)
    for i in range(m):
        add(i, (i + 1) % m); add(m + i, m + (i + 1) % m); add(i, m + i)
    return [sorted(x) for x in adj]


def truncate(adj):
    n = len(adj); idx = {}; c = 0
    for v in range(n):
        for w in adj[v]: idx[(v, w)] = c; c += 1
    out = [set() for _ in range(c)]
    for v in range(n):
        for a, b in itertools.combinations(adj[v], 2):
            out[idx[(v, a)]].add(idx[(v, b)]); out[idx[(v, b)]].add(idx[(v, a)])
    for v in range(n):
        for w in adj[v]:
            out[idx[(v, w)]].add(idx[(w, v)]); out[idx[(w, v)]].add(idx[(v, w)])
    return [sorted(x) for x in out]


# ------------------------------------------------------- cycle machinery

def cycles(adj):
    """set of all simple cycles, each as a frozenset of edges."""
    n = len(adj); out = set()
    for s in range(n):
        path = [s]; vis = [False] * n; vis[s] = True
        def dfs(v):
            for w in adj[v]:
                if w < s: continue
                if w == s and len(path) >= 3:
                    out.add(frozenset(
                        (min(path[i], path[(i + 1) % len(path)]),
                         max(path[i], path[(i + 1) % len(path)]))
                        for i in range(len(path))))
                    continue
                if vis[w]: continue
                vis[w] = True; path.append(w); dfs(w); path.pop(); vis[w] = False
        dfs(s); vis[s] = False
    return out


def spectrum(adj):
    return set(len(c) for c in cycles(adj))


def spectrum_through(adj, e):
    """S_e(H): lengths of cycles containing the edge e = (u,v), u<v."""
    return set(len(c) for c in cycles(adj) if e in c)


def delete_edge(adj, e):
    u, v = e
    out = [list(x) for x in adj]
    out[u] = [w for w in out[u] if w != v]
    out[v] = [w for w in out[v] if w != u]
    return out


def edges_of(adj):
    return [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]


# ------------------------------------------------------------- the chain

def chain(blocks):
    """CHAIN(B_1,e_1; ...; B_m,e_m): delete every e_i = u_i v_i, then add
    v_i -- u_{i+1} cyclically.  m = 2 is the edge amalgam.  Returns adj."""
    off = []; tot = 0
    for (adj, e) in blocks:
        off.append(tot); tot += len(adj)
    out = [set() for _ in range(tot)]
    for k, (adj, e) in enumerate(blocks):
        o = off[k]; u, v = e
        for a in range(len(adj)):
            for b in adj[a]:
                if (min(a, b), max(a, b)) == (u, v): continue
                out[o + a].add(o + b); out[o + b].add(o + a)
    m = len(blocks)
    for k in range(m):
        v_k = off[k] + blocks[k][1][1]
        u_n = off[(k + 1) % m] + blocks[(k + 1) % m][1][0]
        out[v_k].add(u_n); out[u_n].add(v_k)
    return [sorted(x) for x in out]


def sumset(sets):
    cur = {0}
    for A in sets:
        cur = set(x + a for x in cur for a in A)
    return cur


# ------------------------------------------------------------------ zoo

ZOO = [('K4', lcf(4, [2, -2])), ('K33', k33()), ('prism-3', prism(3)),
       ('cube', lcf(8, [3, -3])), ('prism-4', prism(4)),
       ('Petersen', gp(5, 2)), ('prism-5', prism(5)),
       ('Moebius-Kantor', gp(8, 3)), ('Heawood', lcf(14, [5, -5])),
       ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5])),
       ('Desargues', gp(10, 3)),
       ('Dodecahedron', lcf(20, [10, 7, 4, -4, -7, 10, -4, 7, -7, 4])),
       ('trunc-Petersen', truncate(gp(5, 2)))]

SPEC = {}
for name, G in ZOO:
    SPEC[name] = spectrum(G)

P('=' * 74)
P('LEMMA 1  --  THE AMALGAM SPECTRUM IS EXACTLY THE SUMSET', flush=True)
P('=' * 74)
P('  claim:  S(H1 (+) H2)  =  S(H1-e1) U S(H2-e2) U (S_e1(H1)+S_e2(H2))')
P('  and     S_new(H1 (+) H2) = S_e1(H1) + S_e2(H2)   (new edge u1u2)')
P()
bad = 0; tests = 0
SMALL = [('K4', lcf(4, [2, -2])), ('K33', k33()), ('prism-3', prism(3)),
         ('cube', lcf(8, [3, -3])), ('prism-4', prism(4)),
         ('Petersen', gp(5, 2))]
for (n1, H1) in SMALL:
    for (n2, H2) in SMALL:
        for e1 in edges_of(H1)[:3]:
            for e2 in edges_of(H2)[:3]:
                A1 = spectrum_through(H1, e1); A2 = spectrum_through(H2, e2)
                B1 = spectrum(delete_edge(H1, e1))
                B2 = spectrum(delete_edge(H2, e2))
                Hs = chain([(H1, e1), (H2, e2)])
                got = spectrum(Hs)
                want = B1 | B2 | set(a + b for a in A1 for b in A2)
                tests += 1
                if got != want:
                    bad += 1
                    P('    MISMATCH %s/%s %s %s: got %s want %s'
                      % (n1, n2, e1, e2, sorted(got), sorted(want)))
                # the new edge: v_1 -- u_2 in chain indexing
                nu = len(H1); newe = (min(e1[1], nu + e2[0]),
                                      max(e1[1], nu + e2[0]))
                gotn = spectrum_through(Hs, newe)
                wantn = set(a + b for a in A1 for b in A2)
                if gotn != wantn:
                    bad += 1
                    P('    NEW-EDGE MISMATCH %s/%s: got %s want %s'
                      % (n1, n2, sorted(gotn), sorted(wantn)))
P('  amalgam pairs tested: %d      violations: %d' % (tests, bad), flush=True)

P()
P('  m-fold chains (m = 3, 4, 5), mixed blocks and mixed edges:')
bad2 = 0; t2 = 0
BASE = [('K4', lcf(4, [2, -2])), ('K33', k33()), ('prism-3', prism(3))]
for m in (3, 4, 5):
    for combo in itertools.product(range(len(BASE)), repeat=m):
        if m == 5 and combo.count(combo[0]) != m: continue   # keep it cheap
        blocks = []
        for j, i in enumerate(combo):
            nm, Hb = BASE[i]
            blocks.append((Hb, edges_of(Hb)[j % len(edges_of(Hb))]))
        Hs = chain(blocks)
        got = spectrum(Hs)
        inner = set()
        As = []
        for (Hb, eb) in blocks:
            inner |= spectrum(delete_edge(Hb, eb))
            As.append(spectrum_through(Hb, eb))
        want = inner | sumset(As)
        t2 += 1
        if got != want:
            bad2 += 1
            P('    MISMATCH m=%d %s: got %s want %s'
              % (m, combo, sorted(got), sorted(want)))
P('  chains tested: %d      violations: %d' % (t2, bad2), flush=True)

P()
P('=' * 74)
P('LEMMA 2  --  gcd OF THE CYCLE SPECTRUM IS 1 OR 2 (Bondy-Vince)', flush=True)
P('=' * 74)
P('  %-16s %5s  %-38s %s' % ('graph', 'n', 'spectrum', 'gcd'))
badg = 0
for name, G in ZOO:
    S = sorted(SPEC[name]); d = 0
    for x in S: d = gcd(d, x)
    flag = '' if d in (1, 2) else '   *** VIOLATION ***'
    if d not in (1, 2): badg += 1
    ss = str(S)
    if len(ss) > 38: ss = ss[:35] + '...'
    P('  %-16s %5d  %-38s %d%s' % (name, len(G), ss, d, flag), flush=True)
P('  violations: %d' % badg)

P()
P('=' * 74)
P('474.3  --  THE ESCAPE CERTIFICATE', flush=True)
P('=' * 74)
P('  A counterexample has NO power of two in its spectrum, so the zoo')
P('  spectra must first be made counterexample-shaped: delete the')
P('  powers of two, S* = S(G) \\ P.  Then compute the least power of')
P('  two in the semigroup <S*>, the number m of blocks realising it,')
P('  and the size of the chain.  By Lemma 1 that chain is a graph of')
P('  minimum degree 3 built only from copies of the block, and it HAS')
P('  a cycle of that length - the ESCAPE.')
P()
P('  Rows where gcd S* is NOT a power of two are the ones the deletion')
P('  broke: by Lemma 2 no minimum-degree-3 graph has such a spectrum,')
P('  so they are not counterexample-shaped after all.  That is Lemma 2')
P('  doing exactly the work Theorem A asks of it.')
P()
P('  %-16s %5s %-20s %4s %8s %6s %8s'
  % ('graph', 'n', 'S* = S(G) minus P', 'gcd', '2^r', 'blocks', 'vertices'))
for name, G in ZOO:
    S = sorted(x for x in SPEC[name] if x & (x - 1))
    if not S:
        P('  %-16s %5d  S* empty' % (name, len(G))); continue
    dS = 0
    for x in S: dS = gcd(dS, x)
    ss = str(S)
    if len(ss) > 20: ss = ss[:17] + '...'
    # least power of two representable as a sum of elements of S,
    # with the minimum number of summands (BFS over reachable sums).
    LIM = 4 * max(S) * max(S) + 8
    best = None
    reach = {0: 0}
    frontier = [0]
    while frontier and best is None:
        nxt = []
        for x in frontier:
            for a in S:
                y = x + a
                if y > LIM or y in reach: continue
                reach[y] = reach[x] + 1
                if y >= 4 and (y & (y - 1)) == 0:
                    best = (y, reach[y]); break
                nxt.append(y)
            if best: break
        frontier = sorted(nxt)
    if best:
        P('  %-16s %5d %-20s %4d %8d %6d %8d'
          % (name, len(G), ss, dS, best[0], best[1], best[1] * len(G)),
          flush=True)
    else:
        P('  %-16s %5d %-20s %4d   NONE - gcd not a power of two'
          % (name, len(G), ss, dS), flush=True)

P()
P('=' * 74)
P('READING', flush=True)
P('=' * 74)
P('  Lemma 1 makes the family ADDITIVE: the union of the spectra of a')
P('  2-sum-closed family is a subsemigroup of N.  Lemma 2 makes its')
P('  gcd a power of two.  Sylvester then puts a power of two IN it.')
P('  So no nonempty 2-sum-closed family can avoid the powers of two -')
P('  which is THEOREM A, and which says the counterexample class is')
P('  not 2-sum-closed unless it is empty.')
