"""canon401.py -- feasibility audit of the canonical form's numeric skeleton (400.3).

Variables (one-visit straddle template, opposed-pair world):
  A = |O_1|, B = |O_2| odd >= 5;  coordinates on O_1 - {b_1}: u = 0 .. u' = A-2.
  corridor feet s1 < s2 (transitions), same parity (400.2);
  window attachments t1 < s1  and  t2 > s2 (straddle, 400.2);
  runs r1, r2 >= 1 (w->y1, y2->w'); transition lengths m1, m2 >= 2 (or 1 = direct);
  constants c12, c13, c23 in {0, 2}.
Proved-necessary constraints imposed:
  P1  |S_3| = 4 + r1 + r2 + (t2 - t1) odd and >= 5             (392.1, girth)
  P2  transitions: m_i == (A+B)/2 (mod 2)                       (395.3)
  P3  4-combo pin per transition: (1+s_i) + (1+q_i) == c12 - m_i (mod 4), where the
      O_2-foot position q_i (from x, 0..B-2) has its even-completing arc 1+q_i; the
      even combo exists for exactly the parity-matched arc pair - encoded as: there
      is an O_2 position q_i in [0, B-2] with 1+s_i + m_i + 1+q_i even and the pin.
  P4  access pins: 1 + t1 == c13 - r1 - 2  and  1 + (A-2-t2) == c13 - r2 - 2 (mod 4)
                                                                (397.1)
  P5  cross pins {2,3}: for each i: exists arc choice on O_2 with
      (1+q_i) + m_i + (s_i - t1) + 1 + (r1 + 1) == c23 (mod 4)  (392.3-form)
  P6  distances: d(b1,b3) >= 5: (1 + t1) + 1 + r1 + 1 >= 5 and
      (1 + (A-2-t2)) + 1 + r2 + 1 >= 5;  d(b1,b2) >= 5: 1 + s1 + m1 + q1-ish + 1 >= 5
      conservative: 1 + s1 + m1 + 1 + 1 >= 5;  d(b2,b3) >= 5: m1 + s1 - t1 + r1 >= 2
      conservative slack.
  P7  sigma law: c12 + c13 + c23 == 2*sigma (mod 4) for some sigma in {0,1}  (325.3)
Search all assignments in small ranges; report feasible/infeasible per (A,B).
"""
import itertools, sys

feas={}
for A in (5,7,9,11):
    for B in (5,7,9,11):
        pistar=((A+B)//2)%2
        found=None
        for c12,c13,c23 in itertools.product((0,2),repeat=3):
            if (c12+c13+c23)%4 not in (0,2): continue   # P7: 2*sigma in {0,2}
            for t1 in range(0,A-1):
                for s1 in range(t1+1,A-1):
                    for s2 in range(s1,A-1):
                        if (s2-s1)%2: continue          # feet one parity
                        for t2 in range(s2+1,A-1):
                            lam=t2-t1
                            for r1 in range(1,13):
                                # P4 access pins
                                if (1+t1)%4 != (c13-r1-2)%4: continue
                                for r2 in range(1,13):
                                    if (1+(A-2-t2))%4 != (c13-r2-2)%4: continue
                                    S3=4+r1+r2+lam
                                    if S3%2==0 or S3<5: continue    # P1
                                    # P6 distances
                                    if (1+t1)+1+r1+1<5: continue
                                    if (1+(A-2-t2))+1+r2+1<5: continue
                                    ok=False
                                    for m1 in range(1,13):
                                        if m1%2!=pistar: continue   # P2
                                        if 1+s1+m1+2<5: continue    # P6 cons.
                                        for q1 in range(0,B-1):
                                            if ((1+s1)+m1+(1+q1))%2: continue
                                            if ((1+s1)+(1+q1))%4 != (c12-m1)%4: continue  # P3
                                            # P5 cross pin, either arc on O_2
                                            v=( (1+q1)+m1+(s1-t1)+1+(r1+1) )%4
                                            v2=( (1+(B-2-q1))+m1+(s1-t1)+1+(r1+1) )%4
                                            if v==c23%4 or v2==c23%4:
                                                ok=True; break
                                        if ok: break
                                    if ok and found is None:
                                        found=(c12,c13,c23,t1,s1,s2,t2,r1,r2,S3)
        feas[(A,B)]=found
        print('A=%2d B=%2d  %s'%(A,B,'FEASIBLE '+str(found) if found else 'INFEASIBLE'))
