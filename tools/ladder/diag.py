# -*- coding: utf-8 -*-
"""Diagnostics for a PROOF of (PSX): two port-straddling outside chords
always cross.  Measures the mechanisms that could deliver it."""

import sys
import random
from itertools import combinations
from collections import Counter
from nlutil import (ov, enum_matchings, closed_adj, girth, residual,
                    planar_hamiltonian)


def zones(T, g, al, de, be):
    z = [0, 0, 0, 0, 0]
    for t in T:
        if t <= g:
            z[0] += 1
        elif t < al:
            z[1] += 1
        elif t < de:
            z[2] += 1
        elif t < be:
            z[3] += 1
        else:
            z[4] += 1
    return z


def gen_pair(T, g, de, al, be):
    z = zones(T, g, al, de, be)
    for zw in (0, 1):
        for zx in (1, 2):
            for zy in (2, 3):
                for zz in (3, 4):
                    need = Counter([zw, zx, zy, zz])
                    if all(z[k] >= need[k] for k in need):
                        return True
    return False


def gen_any(T, outch):
    for c1, c2 in combinations(outch, 2):
        for (g, de), (al, be) in ((c1, c2), (c2, c1)):
            if g < al < de < be and gen_pair(T, g, de, al, be):
                return True
    return False


def analyse(n, lads):
    st = Counter()
    nopincer = []
    for ch in lads:
        if girth(closed_adj(ch, n)) < 5:
            continue
        st['ladders'] += 1
        pl = planar_hamiltonian(ch)
        if pl:
            st['ladders-PLANAR'] += 1
        for (a, b) in ch:
            R = residual(ch, n, a, b)
            if R is None:
                continue
            m, f, X, ports, outch = R
            st['residual'] += 1
            if pl:
                st['residual-in-PLANAR-ladder'] += 1
            T = sorted(ports)
            F = sorted(f.values())
            muH = {}
            for x, y in outch:
                muH[x] = y
                muH[y] = x
            assert set(muH) | set(F) == set(range(m + 1)), (ch, a, b)
            # 1. pincered F-point
            pin = [t for t in F
                   if t - 1 in muH and t + 1 in muH
                   and muH[t - 1] > t and muH[t + 1] < t]
            if pin:
                st['has-pincer'] += 1
            else:
                st['NO-pincer'] += 1
                nopincer.append((ch, (a, b), m, outch, T, F))
            # 2. straddle counts
            kap = [sum(1 for c in outch if c[0] < t < c[1]) for t in F]
            st['minkappa=%d' % min(kap)] += 1
            # 3. crossing straddlers, and whether each gives GEN
            S = [c for c in outch if any(c[0] < t < c[1] for t in T)]
            xs = [(c, d) for c, d in combinations(S, 2) if ov(c, d)]
            if xs:
                st['straddlers-cross'] += 1
                ok = 0
                for c, d in xs:
                    (g, de), (al, be) = (c, d) if c[0] < d[0] else (d, c)
                    if gen_pair(T, g, de, al, be):
                        ok += 1
                if ok == 0:
                    st['ALL-straddler-crossings-degenerate!'] += 1
                elif ok < len(xs):
                    st['some-straddler-crossing-degenerate'] += 1
            else:
                st['NO-straddler-crossing!'] += 1
            if not gen_any(T, outch):
                st['GEN-FAIL!'] += 1
    return st, nopincer


def main():
    sys.setrecursionlimit(40000)
    rng = random.Random(31415)
    for n, cap in [(14, None), (16, None), (18, 6000), (20, 2500),
                   (22, 1200), (24, 600)]:
        lads = enum_matchings(n, cap, None if cap is None else rng)
        st, nop = analyse(n, lads)
        print("=" * 62)
        print("2N =", n, " matchings:", len(lads))
        for k in sorted(st):
            print("   %-36s %d" % (k, st[k]))
        if nop:
            print("   --- residual instances with NO pincered port:",
                  len(nop))
            for e in nop[:3]:
                print("      lad", e[0], "Z", e[1], "m", e[2])
                print("        outch", e[3], " T", e[4], " F", e[5])


main()
