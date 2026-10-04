import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, sys
from collections import deque
base=_HERE+""
exec(open(base+"surv.py").read().split("def main")[0])

def has2mod4(adj,n,cap=10**7):
    """true if some cycle has length = 2 mod 4 (all cycles, chords allowed)"""
    cnt=0
    for s in range(n):
        stack=[(s,1<<s,0)]
        while stack:
            u,mask,d=stack.pop()
            cnt+=1
            if cnt>cap: return None
            for w in adj[u]:
                if w==s:
                    if d+1>=3 and (d+1)%4==2: return True
                elif w>s and not (mask>>w & 1):
                    stack.append((w,mask|1<<w,d+1))
    return False

def lcf(n,seq):
    E=[(i,(i+1)%n) for i in range(n)]
    S=set((min(a,b),max(a,b)) for a,b in E)
    for i,s in enumerate(seq):
        a,b=i,(i+s)%n
        k=(min(a,b),max(a,b))
        if k not in S: S.add(k); E.append(k)
    adj=[[] for _ in range(n)]
    for a,b in S: adj[a].append(b); adj[b].append(a)
    return adj

cages=[("Petersen",10,[5]*10),
       ("Heawood",14,[5,-5]*7),
       ("Mobius-Kantor",16,[5,-5]*8),
       ("Pappus",18,[5,7,-7,7,-7,-5]*3),
       ("Desargues",20,[5,-5,9,-9]*5),
       ("Dodecahedron",20,[10,7,4,-4,-7,10,-4,7,-7,4]*2),
       ("McGee",24,[12,7,-7]*8),
       ("Nauru",24,[5,-9,7,-7,9,-5]*4),
       ("Tutte-Coxeter",30,[-13,-9,7,-7,9,13]*5)]
print("A) named cubic cages / symmetric graphs")
for nm,n,seq in cages:
    adj=lcf(n,seq)
    g=girth(adj,n); r=has2mod4(adj,n)
    print("   %-15s n=%2d girth=%d   cycle = 2 mod 4 ? %s"%(nm,n,g,r))

print("\nB) small cubic graphs of girth 3 or 4 (known failures)")
for nm,n,E in [("K4",4,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]),
               ("prism",6,[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)]),
               ("K33",6,[(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)]),
               ("cube",8,[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)])]:
    adj=[[] for _ in range(n)]
    for a,b in E: adj[a].append(b); adj[b].append(a)
    print("   %-6s n=%d girth=%d   cycle = 2 mod 4 ? %s"%(nm,n,girth(adj,n),has2mod4(adj,n)))

print("\nC) random cubic, girth >= 5, hunting for a counterexample")
rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 13)
tot=0; bad=0
for n in (10,12,14,16,18,20,22,24):
    for _ in range(80):
        g=randcubic(n,rng)
        if not g: continue
        E,adj=g
        r=has2mod4(adj,n)
        if r is None: continue
        tot+=1
        if not r: bad+=1; print("   COUNTEREXAMPLE n=%d E=%s"%(n,E))
print("   cubic girth>=5 tested:",tot,"  without a 2-mod-4 cycle:",bad)
