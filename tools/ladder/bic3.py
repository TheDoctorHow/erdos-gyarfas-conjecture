import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random
base=_HERE+""
exec(open(base+"bicycle.py").read().split("def main")[0])
src=open(base+"settle.py").read()
exec(src[src.index("def named"):src.index("def build")])
exec(src[src.index("def build"):src.index("def cycles")])
exec(src[src.index("def randcubic"):src.index("def build")])

def bic(n,E):
    E=sorted(set((min(a,b),max(a,b)) for a,b in E))
    cut,cyc=spaces(E,n); mu=len(E)-n+1
    return mu,dim_inter(cut,cyc)

print("named cubic graphs (girth>=5 unless noted)")
for nm,(n,E) in named().items():
    mu,d=bic(n,E)
    print("  %-14s n=%2d mu=%2d  d=%d   d>=mu-2? %s"%(nm,n,mu,d,d>=mu-2))

print("\nrandom CUBIC girth>=5 graphs")
rng=random.Random(5); mx=0
for n in (14,16,18,20,22,24,26,28,30):
    best=0; cnt=0
    for k in range(40):
        E=randcubic(n,rng)
        if not E: continue
        mu,d=bic(n,E); best=max(best,d); cnt+=1; mx=max(mx,d)
    print("  n=%2d  samples=%2d  mu=%2d  max d=%d"%(n,cnt,n//2+1,best))
print("\noverall max bicycle dimension seen:",mx)
