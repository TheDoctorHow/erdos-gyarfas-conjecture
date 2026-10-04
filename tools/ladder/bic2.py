import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random
from collections import deque
exec(open(_HERE+"bicycle.py").read().split("def main")[0])

def main2():
    rng=random.Random(99)
    print("subcubic 2-connected girth>=5, >=1 interior vertex")
    print(" mu   count   max_d   #(d>=mu-2)")
    from collections import defaultdict
    stat=defaultdict(lambda:[0,0,0])
    for n in range(8,27):
        for s in range(3,min(n,10)+1):
            for _ in range(1500):
                g=gen(n,s,rng)
                if not g: continue
                E,adj=g
                if n-s<1: continue
                mu=len(E)-n+1
                cut,cyc=spaces(E,n)
                d=dim_inter(cut,cyc)
                st=stat[mu]; st[0]+=1; st[1]=max(st[1],d)
                if d>=mu-2: st[2]+=1
    for mu in sorted(stat):
        c,md,v=stat[mu]
        print("%3d %7d %7d %10d"%(mu,c,md,v))

    # named cubic graphs
    print()
    for nm,(n,E) in named().items():
        adj=build(n,E)
        cut,cyc=spaces(sorted((min(a,b),max(a,b)) for a,b in E),n)
        mu=len(E)-n+1
        d=dim_inter(cut,cyc)
        print("  %-14s n=%2d mu=%2d  bicycle dim d=%d   d>=mu-2? %s"%(nm,n,mu,d,d>=mu-2))

main2()
