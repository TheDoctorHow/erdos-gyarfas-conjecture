import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools, sys
from collections import deque, Counter
base=_HERE+""
src=open(base+"interrogate.py").read()
exec(src[:src.rindex("def main")])

def main():
    rng=random.Random(314)
    tested=0; single=0; formula=0; bad=[]
    for nH in (18,20):
        for trial in range(160):
            g=randcubic(nH,rng)
            if not g: continue
            E,adj0=g
            D=[dist(adj0,nH,v) for v in range(nH)]
            def edist(e,f): return min(D[a][b] for a in e for b in f)
            cand=[t for t in itertools.combinations(E,3)
                  if all(edist(t[i],t[j])>=3 for i in range(3) for j in range(i+1,3))]
            if not cand: continue
            rng.shuffle(cand)
            for M in cand[:2]:
                adj,term=build_B(nH,E,M); N=nH+len(M)
                for b1,b2 in itertools.permutations(term,2):
                    if two_disj_even(adj,N,b1,b2): continue
                    a=shortest_even(adj,N,b1,b2)
                    if a is None: continue
                    L=len(a)-1; inter=set(a[1:-1])
                    Bpp=set(range(N))-inter-{b2}          # THE REPAIR: delete b2 too
                    if b1 not in Bpp: continue
                    pos={}
                    for i in range(1,L):
                        for u in adj[a[i]]:
                            if u not in a and u in Bpp: pos.setdefault(u,i)
                    for u,i in sorted(pos.items()):
                        S=all_path_lens(adj,Bpp,b1,u)
                        if S is None: continue
                        tested+=1
                        ev={v%4 for v in S if v%2==0}
                        od={v%4 for v in S if v%2==1}
                        cls = ev if ev else od
                        if len(cls)==1: single+=1
                        else: bad.append((L,i,sorted(ev),sorted(od)))
                        # formula: |Q| = i-1 mod 4, and even iff i odd
                        pred=(i-1)%4
                        okpar = (bool(ev) and i%2==1) or (bool(od) and i%2==0)
                        if len(cls)==1 and list(cls)[0]==pred and okpar: formula+=1
    print("b1-u_i pairs tested inside B'' = B - int(alpha) - b2 :",tested)
    print("  with a SINGLE class mod 4                          :",single)
    print("  matching the Position Formula |Q| = i-1 mod 4      :",formula)
    print("  violations                                         :",len(bad))
    for b in bad[:8]: print("    L=%d i=%d even=%s odd=%s"%b)

main()
