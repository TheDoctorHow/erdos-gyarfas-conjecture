"""w294.py -- verify and upgrade the counterexample to 293.4.

293.4 claims P = alpha[b_1,x] + E[x,y] + alpha[y,b_2] is a path, where x is the FIRST
and y the LAST vertex of an odd-shift ear E of beta on the even arc alpha.  It need not
be simple.  This script (1) checks the witness, (2) completes it to a graph in which
every vertex has degree 3 except three terminals of degree 2, girth >= 5, 2-connected,
non-bipartite -- i.e. a graph shaped exactly like a bridge B.
"""
import itertools, sys
from collections import deque

# ---------------------------------------------------------------- witness
alpha = ['b1'] + ['a%d' % i for i in range(1, 10)] + ['b2']          # length 10, EVEN
beta  = ['b1'] + ['c%d' % i for i in range(1, 7)] + ['b2']           # length  7, ODD
E     = ['c2', 'd1', 'd2', 'a5', 'a4', 'e1', 'e2', 'e3',
         'a8', 'a9', 'g1', 'g2', 'g3', 'c6']                         # ear of beta

def pedges(p):
    return [frozenset((p[i], p[i + 1])) for i in range(len(p) - 1)]

edges = set(pedges(alpha)) | set(pedges(beta)) | set(pedges(E))
verts = sorted({v for e in edges for v in e})

def adj_of(es):
    a = {v: set() for v in {x for e in es for x in e}}
    for e in es:
        u, v = tuple(e)
        a[u].add(v); a[v].add(u)
    return a

def girth(a):
    best = 10 ** 9
    for s in a:
        dist = {s: 0}; par = {s: None}; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1; par[w] = u; q.append(w)
                elif w != par[u]:
                    best = min(best, dist[u] + dist[w] + 1)
    return best

def connected(a):
    s = next(iter(a)); seen = {s}; q = deque([s])
    while q:
        u = q.popleft()
        for w in a[u]:
            if w not in seen:
                seen.add(w); q.append(w)
    return len(seen) == len(a)

def two_connected(a):
    if len(a) < 3 or not connected(a):
        return False
    for v in a:
        sub = {u: (a[u] - {v}) for u in a if u != v}
        if not connected(sub):
            return False
    return True

def bipartite(a):
    col = {}
    for s in a:
        if s in col:
            continue
        col[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in col:
                    col[w] = 1 - col[u]; q.append(w)
                elif col[w] == col[u]:
                    return False
    return True

def dist(a, s, t):
    d = {s: 0}; q = deque([s])
    while q:
        u = q.popleft()
        if u == t:
            return d[u]
        for w in a[u]:
            if w not in d:
                d[w] = d[u] + 1; q.append(w)
    return None

A = adj_of(edges)
print('--- witness as given ---')
print('vertices %d  edges %d  maxdeg %d  girth %d' %
      (len(A), len(edges), max(len(A[v]) for v in A), girth(A)))
Q = beta[:beta.index('c2') + 1] + E[1:-1] + beta[beta.index('c6'):]
print('Q simple:', len(Q) == len(set(Q)), ' |Q| =', len(Q) - 1,
      '(EVEN)' if (len(Q) - 1) % 2 == 0 else '(ODD)')
span = beta.index('c6') - beta.index('c2')
print('ear shift sigma = %d - %d = %d  (%s)' %
      (len(E) - 1, span, len(E) - 1 - span, 'ODD' if (len(E) - 1 - span) % 2 else 'EVEN'))
cross = [v for v in E if v in alpha]
print('E n alpha (E-order):', cross, ' -> x =', cross[0], ' y =', cross[-1])
q_x, q_y = alpha.index(cross[0]), alpha.index(cross[-1])
print('q_x = %d, q_y = %d  -> %s' % (q_x, q_y, 'PARALLEL' if q_x < q_y else 'ANTI-PARALLEL'))
m_xy = E.index(cross[-1]) - E.index(cross[0])
print('m_xy = %d ; condition (iii)  q_x-q_y == m_xy (mod 2): %s'
      % (m_xy, (q_x - q_y - m_xy) % 2 == 0))
P293 = alpha[:q_x + 1] + E[E.index(cross[0]) + 1:E.index(cross[-1]) + 1] + alpha[q_y + 1:]
print('293.4 P =', ' '.join(P293))
print('P SIMPLE?', len(P293) == len(set(P293)),
      ' repeated:', sorted({v for v in P293 if P293.count(v) > 1}))

# ---------------------------------------------------------------- upgrade
# Every non-terminal must reach degree 3.  b1, b2 stay at degree 2; choose one further
# degree-2 vertex as the third terminal b3.  Add a matching on the rest.
need = [v for v in verts if len(A[v]) == 2 and v not in ('b1', 'b2')]
print('\n--- upgrade: %d vertices need a third edge ---' % len(need))

import random
best = None
A0 = adj_of(edges)
cands = [v for v in need if dist(A0, 'b1', v) >= 5 and dist(A0, 'b2', v) >= 5]
print('b3 candidates at base distance >= 5 from both b1,b2:', cands)
for b3 in cands or need:
    pool = [v for v in need if v != b3]
    if len(pool) % 2:
        continue
    for trial in range(400):
        rng = random.Random(1000 * cands.index(b3) if b3 in cands else 0)
        rng.seed(trial)
        cur = set(edges)
        used = set()
        order = pool[:]

        def ok(u, v, es):
            if u == v or frozenset((u, v)) in es:
                return False
            return girth(adj_of(es | {frozenset((u, v))})) >= 5

        def bt(i):
            while i < len(order) and order[i] in used:
                i += 1
            if i == len(order):
                return True
            u = order[i]
            opts = [w for w in order[i + 1:] if w not in used]
            rng.shuffle(opts)
            for v in opts:
                if ok(u, v, cur):
                    cur.add(frozenset((u, v))); used.add(u); used.add(v)
                    if bt(i + 1):
                        return True
                    cur.discard(frozenset((u, v))); used.discard(u); used.discard(v)
            return False

        if not bt(0):
            continue
        a = adj_of(cur)
        if not (two_connected(a) and not bipartite(a) and girth(a) >= 5):
            continue
        d12, d13, d23 = dist(a, 'b1', 'b2'), dist(a, 'b1', b3), dist(a, 'b2', b3)
        if min(d12, d13, d23) >= 5:
            best = (b3, set(cur), a, (d12, d13, d23))
            break
    if best:
        break

if best is None:
    print('no completion found')
    sys.exit(1)

b3, ES, A2, dists = best
print('third terminal b3 =', b3)
print('n = %d, m = %d' % (len(A2), len(ES)))
print('degrees: ', sorted({v: len(A2[v]) for v in A2}.items(), key=lambda kv: kv[1]))
print('deg-2 vertices:', sorted(v for v in A2 if len(A2[v]) == 2))
print('girth = %d, 2-connected = %s, bipartite = %s'
      % (girth(A2), two_connected(A2), bipartite(A2)))
print('terminal distances d(b1,b2)=%s d(b1,b3)=%s d(b2,b3)=%s' % dists)

# the four objects must survive verbatim in the completed graph
for name, p in (('alpha', alpha), ('beta', beta), ('E', E), ('Q', Q)):
    good = all(frozenset((p[i], p[i + 1])) in ES for i in range(len(p) - 1))
    print('%-6s still a path: %s  length %d  %s' %
          (name, good and len(p) == len(set(p)), len(p) - 1,
           'EVEN' if (len(p) - 1) % 2 == 0 else 'ODD'))
print('alpha n beta =', sorted(set(alpha) & set(beta)))
print('O = alpha+beta is a cycle of length %d (ODD=%s)'
      % (len(alpha) - 1 + len(beta) - 1, (len(alpha) + len(beta)) % 2 == 0))
print('int(E) n beta =', sorted((set(E[1:-1]) & set(beta))))
print('E n alpha      =', sorted(set(E) & set(alpha)))
print('293.4 P simple in the completed graph:', len(P293) == len(set(P293)))

with open('tools/w294_bridge.txt', 'w') as f:
    for e in sorted(ES, key=lambda s: sorted(s)):
        u, v = sorted(e)
        f.write('%s %s\n' % (u, v))
print('\nedge list written to tools/w294_bridge.txt')
