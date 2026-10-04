"""pocket405.py -- mod-4 composite audit of the canonical form at A = 9 (404.3).

Sound: only proved-necessary constraints are imposed; INFEASIBLE => exclusion.
Layout (A = 9, positions 0..7 on O_1 - {b_1}):
  feet s1 < s2 = s1+2k same parity; window L (< s1, phi = alpha), R (> s2, phi =
  beta), alpha != beta, |L|+|R| >= 3, >= 1 each; every corridor-interior wrong-parity
  position pocketed; pockets >= 3 attachments, pocket-phi constant (404.1); all 8
  positions' third edges committed consistently (window, feet/transition, pocket).
Mod-2 laws: phi (403.1ii), pocket-phi (404.1), feet one parity (400.2).
Mod-4 pins imposed:
  P3  direct transition pins at each foot (with O_2 feet q_i distinct, m_i == pi*)
  P4  access pins at each window attachment: arc_even(t) == c13 - r_t - 2 (mod 4),
      arc side forced by parity matching (400.2)
  PK1 pocket-{1,2} composite: b_1 -(u-side)- pi_i -X- pi_j -arc- s2 -tau2- O_2 - b_2,
      valid when pi_i in left region, pi_j in right region (segments disjoint):
      (1+p_i) + 2 + rho_ij + (p_j - s2) + m2 + (O_2 even arc) == c12 (mod 4)
  PK2 pocket-{1,3} composite: b_1 -(u-side)- pi_i -X- pi_j -arc- r_att -D3- b_3,
      valid when pi_i left, pi_j right, r_att > p_j (arc p_j..r_att disjoint):
      (1+p_i) + 2 + rho_ij + (r_att - p_j) + 1 + (r_(r_att) + 1) == c13 (mod 4)
  rho_ij free in Z_4 with parity p_i + p_j (404.1); r_t free with parity
  W + phi(t) + p_t.
Search: all layouts, all bits, all mod-4 assignments in range. Report per B.
"""
import itertools, sys

A=9; POS=list(range(A-1))   # 0..7
results={}
for B in (5,7,9,11,13):
    pistar=((A+B)//2)%2
    feasible=None
    for s1 in range(1,A-2):
        for gap in (2,4):
            s2=s1+gap
            if s2>A-3: continue
            feet={s1,s2}
            interior=[p for p in range(s1+1,s2) if (p-s1)%2==1]  # wrong parity
            left=[p for p in range(0,s1)]; right=[p for p in range(s2+1,A-1)]
            # window: choose L subset of left, R subset of right, |L|+|R|>=3
            for Lsz in (1,2,3):
                for Rsz in (1,2,3):
                    if Lsz+Rsz<3: continue
                    for L in itertools.combinations(left,Lsz):
                        for R in itertools.combinations(right,Rsz):
                            used=set(L)|set(R)|feet
                            rest=[p for p in POS if p not in used]
                            # pocket: must contain all interior positions + others from rest
                            if any(p not in rest for p in interior): continue
                            others=[p for p in rest if p not in interior]
                            need=max(0,3-len(interior))
                            for extra in itertools.combinations(others,need):
                                X=sorted(set(interior)|set(extra))
                                spare=[p for p in rest if p not in X]
                                if spare: continue   # every third edge must be committed
                                if len(X)<3: continue
                                # bits and constants
                                ok=False
                                for c12,c13,c23 in itertools.product((0,2),repeat=3):
                                    if (c12+c13+c23)%4 not in (0,2): continue
                                    for W,alpha in itertools.product((0,1),(0,1)):
                                        beta=1-alpha
                                        # access pins: r_t exists for any residue; record forced r_t mod 4
                                        def rt(p,phi):
                                            # arc side: u-side arc 1+p (parity p+1), u'-side A-1-p (parity p)
                                            # arc_even parity must == r parity? even composite: arc == r (mod 2)
                                            rp=(W+phi+p)%2
                                            for arc in (1+p, A-1-p):
                                                if arc%2==rp%2:
                                                    return (c13-arc-2)%4
                                            return None
                                        rmap={}
                                        good=True
                                        for p in L:
                                            v=rt(p,alpha)
                                            if v is None: good=False; break
                                            rmap[p]=v
                                        if good:
                                            for p in R:
                                                v=rt(p,beta)
                                                if v is None: good=False; break
                                                rmap[p]=v
                                        if not good: continue
                                        # transitions at feet: find (m,q) pairs
                                        def tok(si):
                                            outs=[]
                                            for m in range(1,13):
                                                if m%2!=pistar: continue
                                                for q in range(0,B-1):
                                                    if ((1+si)+m+(1+q))%2: continue
                                                    if ((1+si)+(1+q))%4!=(c12-m)%4: continue
                                                    outs.append((m,q))
                                            return outs
                                        o1=tok(s1); o2=tok(s2)
                                        pairs=[(t1,t2) for t1 in o1 for t2 in o2 if t1[1]!=t2[1]]
                                        if not pairs: continue
                                        # pocket composites: choose rho_ij in Z4 with parity p_i+p_j,
                                        # satisfying PK1 (through s2, right-side pi_j) and PK2 (to R attachments)
                                        Xl=[p for p in X if p<s1]; Xr=[p for p in X if p>s2]
                                        Xi=[p for p in X if s1<p<s2]
                                        # constraints only for cross pairs with valid geometry:
                                        # (i in Xl union Xi-with-u-side-access?) -- SOUND subset: i in Xl, j in Xr
                                        cross=[(i,j) for i in Xl for j in Xr]
                                        if not cross:
                                            ok=True; break   # no PK constraints encoded: feasible (sound)
                                        found_rho=False
                                        for rhovec in itertools.product(range(4),repeat=len(cross)):
                                            all_ok=True
                                            for (idx,(i,j)) in enumerate(cross):
                                                rho=rhovec[idx]
                                                if rho%2!=(i+j)%2: all_ok=False; break
                                                # PK1 with transition at s2: need SOME (m2,q2) in o2:
                                                pk1=False
                                                for (m2,q2) in o2:
                                                    base=(1+i)+2+rho+(j-s2)+m2
                                                    for arcB in (1+q2, B-1-q2):
                                                        if (base+arcB)%2==0 and (base+arcB)%4==c12%4:
                                                            pk1=True; break
                                                    if pk1: break
                                                if not pk1: all_ok=False; break
                                                # PK2 to each R attachment beyond j:
                                                for ra in R:
                                                    if ra<=j: continue
                                                    tot=(1+i)+2+rho+(ra-j)+1+(rmap[ra]+1)
                                                    if tot%2==0 and tot%4!=c13%4: all_ok=False; break
                                                    # if odd: not an even connecting path: unconstrained
                                                if not all_ok: break
                                            if all_ok: found_rho=True; break
                                        if found_rho: ok=True; break
                                    if ok: break
                                if ok and feasible is None:
                                    feasible=(s1,s2,L,R,tuple(X))
                            if feasible: break
                        if feasible: break
                    if feasible: break
                if feasible: break
            if feasible: break
        if feasible: break
    results[B]=feasible
    print('A=9 B=%2d  %s'%(B,'FEASIBLE '+str(feasible) if feasible else 'INFEASIBLE'))
