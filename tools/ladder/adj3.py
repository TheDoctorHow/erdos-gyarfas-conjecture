# -*- coding: utf-8 -*-
"""Which (ADJ) witnesses survive adversarial minimisation?  If the
surviving core is always at the same kind of position, that is the
mechanism to prove."""

POW = set(1 << i for i in range(1, 22))

EXTREMAL = {
 20: [(1,5),(2,12),(3,20),(4,13),(6,11),(7,17),(8,14),(9,18),(10,16),
      (15,19)],
 26: [(1,6),(2,12),(3,8),(4,10),(5,17),(7,20),(9,18),(11,15),(13,22),
      (14,25),(16,24),(19,23),(21,26)],
 32: [(1,10),(2,7),(3,11),(4,8),(5,17),(6,23),(9,31),(12,18),(13,26),
      (14,24),(15,20),(16,32),(19,28),(21,30),(22,27),(25,29)],
 40: [(1,15),(2,11),(3,17),(4,13),(5,10),(6,18),(7,37),(8,12),(9,32),
      (14,22),(16,28),(19,27),(20,25),(21,33),(23,34),(24,38),(26,40),
      (29,35),(30,39),(31,36)],
 50: [(1,7),(2,27),(3,8),(4,44),(5,10),(6,12),(9,25),(11,38),(13,29),
      (14,18),(15,24),(16,21),(17,33),(19,23),(20,34),(22,49),(26,31),
      (28,32),(30,35),(36,40),(37,46),(39,45),(41,47),(42,50),(43,48)],
}

def overlap(c, d):
    (a, b), (p, q) = c, d
    return not (b < p or q < a)

def l1(c, d):
    return abs(c[0] - d[0]) + abs(c[1] - d[1])

for n, ch in sorted(EXTREMAL.items()):
    mu = {}
    for a, b in ch:
        mu[a] = b
        mu[b] = a
    print("=== 2N =", n)
    for p in range(1, n):
        c = tuple(sorted((p, mu[p])))
        cp = tuple(sorted((p + 1, mu[p + 1])))
        if c == cp:
            continue
        Q = (min(mu[p], mu[p + 1]), max(mu[p], mu[p + 1]))
        beta = (mu[p] < p) + (mu[p + 1] < p)
        tp = 'L' if mu[p] > p else 'R'
        tq = 'L' if mu[p + 1] > p + 1 else 'R'
        for D in ch:
            if D == c or D == cp:
                continue
            if not overlap(D, Q):
                continue
            if (beta + (D[0] < p) + (D[1] < p)) % 2:
                continue
            w = 4 + l1(Q, D)
            if w in POW:
                print("   p=%2d (%s%s)  Q=%-9s D=%-9s w=%d"
                      % (p, tp, tq, str(Q), str(D), w))
