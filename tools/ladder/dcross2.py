import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
# Trichotomy check (anti-parallel branch) + the parallel branch.
base = _HERE
src = open(base + "dcross.py").read()
exec(src[:src.index("tot = 0")])

# ---- anti-parallel: full trichotomy -------------------------------------
c_a = c_b = c_c = 0
err = 0
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
                                    e1 = r["L1"] % 2 == 0
                                    e2 = r["L2"] % 2 == 0
                                    if r["cond3"] and r["cond4"]:
                                        c_a += 1
                                        if not (e1 and e2):
                                            err += 1
                                    elif r["cond3"]:
                                        c_b += 1
                                        if e1 or e2:      # expect both ODD
                                            err += 1
                                    else:
                                        c_c += 1
                                        if e1 == e2:      # expect exactly one
                                            err += 1
print("anti-parallel branch A  (iii)&(iv)  -> two disjoint EVEN paths :", c_a)
print("anti-parallel branch B  (iii)&~(iv) -> two disjoint ODD  paths :", c_b)
print("anti-parallel branch C  ~(iii)      -> exactly one even        :", c_c)
print("violations of the trichotomy                                  :", err)

# ---- parallel: y AFTER x along alpha ------------------------------------
tot2 = badp = 0
even2 = 0
for A in (4, 6, 8, 10):
    for B in (3, 5, 7, 9):
        for p in range(1, B - 1):
            for pt in range(p + 1, B):
                for qx in range(1, A - 1):
                    for qy in range(qx + 1, A):
                        for m1 in (1, 2, 3):
                            for mxy in (1, 2, 3):
                                for m2 in (1, 2, 3):
                                    a, c, E, s, t, x, y, u, v, w = build(
                                        A, B, p, pt, qy, qx, m1, mxy, m2)
                                    if len(set(E)) != len(E):
                                        continue
                                    if (m1 + mxy + m2 - (pt - p)) % 2 == 0:
                                        continue
                                    tot2 += 1
                                    # shortcut path P = alpha[b1,x]+E[x,y]+alpha[y,b2]
                                    P = a[0:qx + 1] + v + a[qy:A + 1]
                                    if len(set(P)) != len(P):
                                        badp += 1
                                        continue
                                    if set(P) & set(c) != {"b1", "b2"}:
                                        badp += 1
                                        continue
                                    LP = len(P) - 1
                                    cond3 = ((qy - qx) - mxy) % 2 == 0
                                    if (LP % 2 == 0) != cond3:
                                        badp += 1
                                    if cond3:
                                        even2 += 1
print()
print("parallel configurations with odd-shift ear                    :", tot2)
print("  shortcut P a path, internally disjoint from beta, |P| parity")
print("  matching cond (iii) -- violations                           :", badp)
print("  of these, |P| EVEN  -> new odd cycle P+beta through the pair:", even2)
