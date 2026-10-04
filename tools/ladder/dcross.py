# Exhaustive synthetic verification of the DOUBLE-CROSSING LEMMA.
#
# Configuration: odd cycle O through b1,b2 with arcs alpha (even, length A)
# and beta (odd, length B).  An ear E of beta from s=c_p to t=c_pt whose
# interior meets alpha first at x=a_qx and last at y=a_qy.
#
# Claim: if  qy < qx  (anti-parallel),  (qx-qy) == mxy mod 2   (cond iii),
# and  p+m1+qx  even  (cond iv), then
#     P1 = beta[b1,s] + E[s,x] + alpha[x,b2]
#     P2 = alpha[b1,y] + E[y,t] + beta[t,b2]
# are internally disjoint b1-b2 paths, both of EVEN length.

def build(A, B, p, pt, qy, qx, m1, mxy, m2):
    a = ["a%d" % i for i in range(A + 1)]
    c = ["c%d" % i for i in range(B + 1)]
    a[0] = c[0] = "b1"
    a[A] = c[B] = "b2"
    u = ["u%d" % i for i in range(1, m1)]
    v = ["v%d" % i for i in range(1, mxy)]
    w = ["w%d" % i for i in range(1, m2)]
    s, t, x, y = c[p], c[pt], a[qx], a[qy]
    E = [s] + u + [x] + v + [y] + w + [t]
    return a, c, E, s, t, x, y, u, v, w


def check(A, B, p, pt, qy, qx, m1, mxy, m2):
    a, c, E, s, t, x, y, u, v, w = build(A, B, p, pt, qy, qx, m1, mxy, m2)
    # ear must be a path with no repeats, interior off beta
    if len(set(E)) != len(E):
        return None
    if len(set(E[1:-1]) & set(c)) != 0:
        return None
    shift = (m1 + mxy + m2) - (pt - p)
    if shift % 2 == 0:
        return None                       # need ODD shift
    cond3 = ((qx - qy) - mxy) % 2 == 0
    cond4 = (p + m1 + qx) % 2 == 0
    P1 = c[0:p + 1] + u + [x] + a[qx + 1:A + 1]
    P2 = a[0:qy + 1] + w + [t] + c[pt + 1:B + 1]
    L1, L2 = len(P1) - 1, len(P2) - 1
    ok_path = len(set(P1)) == len(P1) and len(set(P2)) == len(P2)
    inter = set(P1) & set(P2)
    ok_disj = inter == {"b1", "b2"}
    return dict(cond3=cond3, cond4=cond4, L1=L1, L2=L2,
                ok_path=ok_path, ok_disj=ok_disj,
                sum_even=(L1 + L2) % 2 == 0)


tot = 0
bad_disj = 0
bad_sum = 0
bad_concl = 0
hits = 0
for A in (4, 6, 8, 10):
    for B in (3, 5, 7, 9):
        for p in range(1, B - 1):
            for pt in range(p + 1, B):
                for qy in range(1, A - 1):
                    for qx in range(qy + 1, A):
                        for m1 in (1, 2, 3):
                            for mxy in (1, 2, 3):
                                for m2 in (1, 2, 3):
                                    r = check(A, B, p, pt, qy, qx,
                                              m1, mxy, m2)
                                    if r is None:
                                        continue
                                    tot += 1
                                    if not (r["ok_path"] and r["ok_disj"]):
                                        bad_disj += 1
                                    # equivalence: sum even  <=>  cond3
                                    if r["sum_even"] != r["cond3"]:
                                        bad_sum += 1
                                    if r["cond3"] and r["cond4"]:
                                        hits += 1
                                        if not (r["L1"] % 2 == 0
                                                and r["L2"] % 2 == 0):
                                            bad_concl += 1

print("configurations with odd-shift ear, anti-parallel crossing :", tot)
print("failures of internal disjointness / simplicity            :", bad_disj)
print("failures of  ( |P1|+|P2| even  <=>  cond (iii) )           :", bad_sum)
print("configurations satisfying (iii) and (iv)                   :", hits)
print("of those, failures of  |P1|,|P2| both even                 :", bad_concl)
