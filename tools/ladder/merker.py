"""Merker's lower-bound family (arXiv:1905.09101, Theorem 7), implemented verbatim.

  H_k : two disjoint cycles  C1 = x_0 x_1 ... x_{2k} x_0  and  C2 = y_0 y_1 ... y_{2k} y_0,
        together with chords  x_i x_{2k-i}  and  y_i y_{2k-i}  for i in {1,...,k-1},
        and one edge  x_k y_k  joining C1 and C2.
        => 4k+2 vertices; the four vertices x_0, x_{2k}, y_0, y_{2k} have degree 2
           (they are the terminals), every other vertex has degree 3.

  D(e,H_k) for e = uv, u ~ u1,u2,v and v ~ v1,v2,u :  delete u and v, add a copy of H_k,
        and add the edges  u1 x_0,  u2 x_{2k},  v1 y_0,  v2 y_{2k}.
        D(M,H_k) applies this successively to every edge of a matching M.

  D_n : cycles C1 = u_1...u_n u_1, C2 = v_1...v_n v_1, C3 = w_1 w_2 ... w_{2n} w_1,
        plus edges  v_i w_{2i-1}  and  u_i w_{2i}  for i in {1,...,n}.   (4n vertices)
        M = { v_i w_{2i-1} } u { u_i w_{2i} }  is a perfect matching (2n edges).

  G(n,k) = D_n(M, H_{k-1}),  n >= 4k+2.
  Theorem 7:  C(G(n,k)) ∩ [2k, 4k+1] = ∅   and   circumference > 2k.

Vertex count: every one of the 4n vertices of D_n is matched, so all 2n matched edges are
replaced; |V(G(n,k))| = 2n * (4(k-1)+2) = 2n(4k-2).
"""
import sys, collections
sys.path.insert(0, '.')
from ref import to_graph6


def build_Dn(n):
    """Return (vertices, edges, matching) for D_n. Vertices are labelled tuples."""
    U = [('u', i) for i in range(1, n + 1)]
    V = [('v', i) for i in range(1, n + 1)]
    W = [('w', j) for j in range(1, 2 * n + 1)]
    E = []
    for i in range(1, n + 1):
        E.append((('u', i), ('u', i % n + 1)))
        E.append((('v', i), ('v', i % n + 1)))
    for j in range(1, 2 * n + 1):
        E.append((('w', j), ('w', j % (2 * n) + 1)))
    M = []
    for i in range(1, n + 1):
        M.append((('v', i), ('w', 2 * i - 1)))
        M.append((('u', i), ('w', 2 * i)))
    return U + V + W, E + M, M


def Hk_edges(k, base):
    """Edges of H_k on vertex ids base..base+4k+1, and its four terminals.
    x_i = base + i (i = 0..2k);  y_i = base + (2k+1) + i."""
    x = lambda i: base + i
    y = lambda i: base + (2 * k + 1) + i
    E = []
    for i in range(2 * k + 1):
        E.append((x(i), x((i + 1) % (2 * k + 1))))
        E.append((y(i), y((i + 1) % (2 * k + 1))))
    for i in range(1, k):
        E.append((x(i), x(2 * k - i)))
        E.append((y(i), y(2 * k - i)))
    E.append((x(k), y(k)))
    return E, (x(0), x(2 * k), y(0), y(2 * k))


def build_G(n, k):
    """G(n,k) = D_n(M, H_{k-1})."""
    j = k - 1                                  # gadget index
    verts, edges, M = build_Dn(n)
    Mset = set(map(lambda e: frozenset(e), M))
    nonM = [e for e in edges if frozenset(e) not in Mset]

    # for each D_n vertex, list its two non-matching incident edges (in a fixed order)
    inc = collections.defaultdict(list)
    for idx, (a, b) in enumerate(nonM):
        inc[a].append(idx)
        inc[b].append(idx)
    for a in verts:
        assert len(inc[a]) == 2, (a, inc[a])

    # one gadget per matched edge; terminals x0,x2k -> u's two non-M edges,
    #                                        y0,y2k -> v's two non-M edges
    port = {}                                   # (nonM edge index, endpoint) -> gadget vertex
    E = []
    nid = 0
    for (u, v) in M:
        ge, (x0, x2k, y0, y2k) = Hk_edges(j, nid)
        nid += 4 * j + 2
        E += ge
        eu = inc[u]; ev = inc[v]
        port[(eu[0], u)] = x0
        port[(eu[1], u)] = x2k
        port[(ev[0], v)] = y0
        port[(ev[1], v)] = y2k
    for idx, (a, b) in enumerate(nonM):
        E.append((port[(idx, a)], port[(idx, b)]))

    deg = collections.Counter()
    for a, b in E:
        deg[a] += 1; deg[b] += 1
    assert nid == 2 * n * (4 * k - 2), (nid, 2 * n * (4 * k - 2))
    assert set(deg.values()) == {3}, sorted(set(deg.values()))
    assert len(set(map(lambda e: tuple(sorted(e)), E))) == len(E), "multi-edge"
    # connectivity
    adj = collections.defaultdict(list)
    for a, b in E:
        adj[a].append(b); adj[b].append(a)
    seen = {0}; st = [0]
    while st:
        z = st.pop()
        for w in adj[z]:
            if w not in seen:
                seen.add(w); st.append(w)
    assert len(seen) == nid, "disconnected"
    return nid, [tuple(sorted(e)) for e in E]


if __name__ == '__main__':
    n = int(sys.argv[1]); k = int(sys.argv[2])
    assert n >= 4 * k + 2, "Theorem 7 requires n >= 4k+2"
    N, E = build_G(n, k)
    print(to_graph6(N, E))
    print("# G(%d,%d): n=%d  m=%d  gadget H_%d  claimed gap [%d,%d]"
          % (n, k, N, len(E), k - 1, 2 * k, 4 * k + 1), file=sys.stderr)
