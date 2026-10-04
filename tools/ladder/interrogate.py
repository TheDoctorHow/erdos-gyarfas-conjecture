import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools, sys
from collections import deque, Counter
base=_HERE+""
src=open(base+"fan.py").read()
exec(src[:src.rindex("def main")])

def all_path_lens(adj,verts,s,t,cap=400000):
    out=set(); stack=[(s,{s},0)]; c=0
    while stack:
        u,seen,d=stack.pop()
        c+=1
        if c>cap: return None
        for w in adj[u]:
            if w not in verts: continue
            if w==t: out.add(d+1)
            elif w not in seen: stack.append((w,seen|{w},d+1))
    return out

def bipclass(adj,verts,root):
    col={root:0}; q=deque([root])
    while q:
        u=q.popleft()
        for w in adj[u]:
            if w not in verts: continue
            if w not in col: col[w]=1-col[u]; q.append(w)
            elif col[w]==col[u]: return None
    return col

def main():
    rng=random.Random(314); shown=0
    for nH in (18,20):
        for trial in range(160):
            if shown>=2: return
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
                if shown>=2: return
                adj,term=build_B(nH,E,M); N=nH+len(M)
                for b1,b2 in itertools.permutations(term,2):
                    if shown>=2: return
                    if two_disj_even(adj,N,b1,b2): continue
                    a=shortest_even(adj,N,b1,b2)
                    if a is None: continue
                    L=len(a)-1; inter=set(a[1:-1]); Bp=set(range(N))-inter
                    T1={b1}; pos={b1:0}
                    for i in range(1,L):
                        if i%2==1:
                            for u in adj[a[i]]:
                                if u not in a: T1.add(u); pos.setdefault(u,i)
                    shown+=1
                    print("="*66)
                    print("INSTANCE %d :  nH=%d  L=%d  alpha=%s"%(shown,nH,L,a))
                    print("  terminals b1=%d b2=%d b3=%d"%(b1,b2,[t for t in term if t not in(b1,b2)][0]))
                    print("  T1 = %s   with alpha-positions %s"%(sorted(T1),{k:pos[k] for k in sorted(T1)}))
                    comp=bipclass(adj,Bp,b1)
                    print("  B' bipartite from b1 ?", "YES" if comp else "NO",
                          " |B'|=%d  connected-from-b1=%d"%(len(Bp),len(comp) if comp else -1))
                    print("\n  SUSPECT (a)+(b): actual path-length classes inside B'")
                    for x,y in itertools.combinations(sorted(T1),2):
                        S=all_path_lens(adj,Bp,x,y)
                        if S is None: print("    %2d-%2d : (too many)"%(x,y)); continue
                        ev=sorted({v%4 for v in S if v%2==0})
                        od=sorted({v%4 for v in S if v%2==1})
                        i,j=pos[x],pos[y]
                        pred = (j-1)%4 if i==0 else (j-i-2)%4
                        print("    %2d-%2d (pos %d,%d): even classes %s  odd classes %s   PREDICTED %d"
                              %(x,y,i,j,ev,od,pred))
                    print("\n  SUSPECT (c): fan-carriers and their bipartition classes")
                    good=[]
                    for w in sorted(Bp):
                        if w in T1 or len(adj[w])!=3: continue
                        if maxflow_to_set(adj,Bp,w,T1,3)>=3: good.append(w)
                    print("    fan-carriers:",good)
                    if comp: print("    their classes:",[comp.get(w,'?') for w in good])
                    adjpairs=[(w,x) for w in good for x in adj[w] if x in good and w<x]
                    print("    adjacent fan-carrying pairs:",adjpairs)

main()
