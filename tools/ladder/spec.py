import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools, sys
from collections import deque, Counter
base=_HERE+""
exec(open(base+"surv.py").read().split("def main")[0])

def allcycles(adj,n):
    """every cycle as a vertex tuple (subcubic: cycles are induced by vertex sets)"""
    out=set()
    for s in range(n):
        stack=[(s,(s,),1<<s)]
        while stack:
            u,pth,mask=stack.pop()
            for w in adj[u]:
                if w==s and len(pth)>=3: out.add(frozenset(pth))
                elif w>s and not (mask>>w & 1):
                    stack.append((w,pth+(w,),mask|1<<w))
    return out

def cyclens(adj,n):
    L=Counter()
    for C in allcycles(adj,n):
        k=len(C)
        # verify 2-regular induced
        if all(sum(1 for w in adj[v] if w in C)==2 for v in C): L[k]+=1
    return L

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 7)
    print("A) cubic graphs of girth>=5: is there a cycle of length = 2 mod 4 ?")
    miss=0; tot=0
    for n in (10,12,14,16,18,20):
        for _ in range(60):
            g=randcubic(n,rng)
            if not g: continue
            E,adj=g
            L=cyclens(adj,n); tot+=1
            if not any(k%4==2 for k in L):
                miss+=1
                print("   NO 2-mod-4 CYCLE: n=%d spectrum=%s"%(n,sorted(L)))
    print("   cubic girth-5 graphs tested:",tot," without a 2-mod-4 cycle:",miss)

    print("\nB) named small cubic graphs")
    for nm,n,E in [("K4",4,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]),
                   ("K33",6,[(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)]),
                   ("Petersen",10,[(0,1),(1,2),(2,3),(3,4),(4,0),(5,7),(7,9),(9,6),(6,8),(8,5),
                                   (0,5),(1,6),(2,7),(3,8),(4,9)])]:
        adj=[[] for _ in range(n)]
        for a,b in E: adj[a].append(b); adj[b].append(a)
        L=cyclens(adj,n)
        print("   %-9s girth=%d  lengths=%s   has 2 mod 4: %s"%(
            nm,min(L),sorted(L),any(k%4==2 for k in L)))

    print("\nC) surviving bridges B (cubic H girth>=5 + far-apart matching, s=3)")
    miss=0; tot=0
    for nH in (18,20):
        for _ in range(60):
            g=randcubic(nH,rng)
            if not g: continue
            E,adj=g
            D=[dist(adj,nH,v) for v in range(nH)]
            def edist(e,f): return min(D[a][b] for a in e for b in f)
            cand=[t for t in itertools.combinations(E,3)
                  if all(edist(t[i],t[j])>=3 for i in range(3) for j in range(i+1,3))]
            if not cand: continue
            M=rng.choice(cand); Ms=set(M)
            L=Counter()
            for C in allcycles(adj,nH):
                ce=[(u,w) for u in C for w in adj[u] if w in C and u<w]
                if len(ce)!=len(C): continue
                L[len(C)+sum(1 for e in ce if e in Ms)]+=1
            tot+=1
            if not any(k%4==2 for k in L):
                miss+=1; print("   NO 2-mod-4 CYCLE in B: nH=%d"%nH)
    print("   surviving bridges tested:",tot," without a 2-mod-4 cycle:",miss)

main()
