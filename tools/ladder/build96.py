"""Build the explicit m=1 K3,3-template candidate (n=96) from the derived design.

Corridors c[x][y] from u_x to v_y, length 11 (10 interiors each).
Interior offset k (u-frame, 1..10) = k-th vertex from u_x along the corridor.
Attachment design: pool template A1B1 A5B5 B2C2 B4C4 A3C5 per branch vertex,
exports {A2,A4,B3,C1,C3}, one mixed edge, 14 long-range export pairings.
"""
import subprocess, sys

# vertex numbering: u1,u2,u3 = 0,1,2 ; v1,v2,v3 = 3,4,5
# corridor interiors: id(x,y,k) for x,y in 1..3, k in 1..10 (u-frame offset)
def iv(x, y, k):
    assert 1 <= k <= 10
    return 6 + ((x - 1) * 3 + (y - 1)) * 10 + (k - 1)

N = 96
edges = set()

def add(a, b):
    e = (min(a, b), max(a, b))
    assert e[0] != e[1], e
    assert e not in edges, ("dup", e)
    edges.add(e)

# skeleton corridors
for x in (1, 2, 3):
    for y in (1, 2, 3):
        u = x - 1
        v = 3 + (y - 1)
        add(u, iv(x, y, 1))
        for k in range(1, 10):
            add(iv(x, y, k), iv(x, y, k + 1))
        add(iv(x, y, 10), v)

# u-frame offset helper for v-frame: v-offset j  ->  u-offset 11 - j
def vf(x, y, j):
    return iv(x, y, 11 - j)

att = []

def A(pq):
    att.append(pq)

# ---- pools: template edges (local frames) ----
# u1: A=c11 B=c12 C=c13   (u-frame = local frame)
A((iv(1,1,1), iv(1,2,1))); A((iv(1,1,5), iv(1,2,5)))
A((iv(1,2,2), iv(1,3,2))); A((iv(1,2,4), iv(1,3,4)))
A((iv(1,1,3), iv(1,3,5)))
# u2: A=c22 B=c23 C=c21
A((iv(2,2,1), iv(2,3,1))); A((iv(2,2,5), iv(2,3,5)))
A((iv(2,3,2), iv(2,1,2))); A((iv(2,3,4), iv(2,1,4)))
A((iv(2,2,3), iv(2,1,5)))
# u3: A=c33 B=c31 C=c32
A((iv(3,3,1), iv(3,1,1))); A((iv(3,3,5), iv(3,1,5)))
A((iv(3,1,2), iv(3,2,2))); A((iv(3,1,4), iv(3,2,4)))
A((iv(3,3,3), iv(3,2,5)))
# v1: A=c11 B=c21 C=c31  (local = v-frame)
A((vf(1,1,1), vf(2,1,1))); A((vf(1,1,5), vf(2,1,5)))
A((vf(2,1,2), vf(3,1,2))); A((vf(2,1,4), vf(3,1,4)))
A((vf(1,1,3), vf(3,1,5)))
# v2: A=c22 B=c32 C=c12
A((vf(2,2,1), vf(3,2,1))); A((vf(2,2,5), vf(3,2,5)))
A((vf(3,2,2), vf(1,2,2))); A((vf(3,2,4), vf(1,2,4)))
A((vf(2,2,3), vf(1,2,5)))
# v3: A=c33 B=c13 C=c23
A((vf(3,3,1), vf(1,3,1))); A((vf(3,3,5), vf(1,3,5)))
A((vf(1,3,2), vf(2,3,2))); A((vf(1,3,4), vf(2,3,4)))
A((vf(3,3,3), vf(2,3,5)))

# ---- mixed edge at v1: (c31@u10, c21@u3) ; v-frames (1, 8) ----
A((iv(3,1,10), iv(2,1,3)))

# ---- long-range odd pairs (u-frame odd offsets) ----
A((iv(1,2,3), iv(2,3,3)))
A((iv(1,3,1), iv(2,1,1)))
A((iv(1,3,3), iv(3,2,1)))
A((iv(3,1,3), iv(2,2,9)))
A((iv(3,2,3), iv(1,1,9)))
A((iv(1,1,7), iv(3,3,9)))
A((iv(2,2,7), iv(3,3,7)))
# ---- long-range even pairs ----
A((iv(1,1,2), iv(2,3,8)))
A((iv(1,1,4), iv(3,2,8)))
A((iv(2,2,2), iv(1,3,8)))
A((iv(2,2,4), iv(3,1,8)))
A((iv(3,3,2), iv(1,2,8)))
A((iv(3,3,4), iv(2,1,8)))
A((iv(1,2,10), iv(2,3,10)))

assert len(att) == 45, len(att)
for p, q in att:
    add(p, q)

# degree check
deg = [0] * N
for a, b in edges:
    deg[a] += 1
    deg[b] += 1
bad = [(i, d) for i, d in enumerate(deg) if d != 3]
print("edges:", len(edges), " non-cubic vertices:", bad)

# graph6 output
adj = [[0] * N for _ in range(N)]
for a, b in edges:
    adj[a][b] = adj[b][a] = 1
bits = []
for j in range(1, N):
    for i in range(j):
        bits.append(adj[i][j])
while len(bits) % 6:
    bits.append(0)
out = [chr(126), chr(((N >> 12) & 63) + 63), chr(((N >> 6) & 63) + 63), chr((N & 63) + 63)]
for k in range(0, len(bits), 6):
    v = 0
    for b in bits[k:k + 6]:
        v = (v << 1) | b
    out.append(chr(v + 63))
g6 = "".join(out)
with open("cand96.g6", "w") as f:
    f.write(g6 + "\n")
print("graph6 written: cand96.g6")
