import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools, sys
from collections import deque, Counter
base=_HERE+""
exec(open(base+"surv.py").read().split("def main")[0])

def build_B(nH,E,M):
    Ms=set(M); adj=[[] for _ in range(nH+len(M))]; term=[]; nxt=nH
    for (a,b) in E:
        if (a,b) in Ms:
            adj[a].append(nxt); adj[nxt].append(a)
            adj[b].append(nxt); adj[nxt].append(b); term.append(nxt); nxt+=1
        else:
            adj[a].append(b); adj[b].append(a)
    return adj,term

def even_paths(adj,N,s,t,cap=200000):
    """all even s-t paths; return list of vertex-sets"""
    out=[]; stack=[(s,1<<s,(s,))]; c=0
    while stack:
        u,mask,p=stack.pop()
        c+=1
        if c>cap: return None
        for w in adj[u]:
            if w==t:
                if (len(p))%2==0: out.append(p+(t,))
            elif not (mask>>w & 1): stack.append((w,mask|1<<w,p+(w,)))
    return out

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 123)
    stat=Counter(); detail=[]
    for nH in (18,20):
        for trial in range(90):
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
                for b1,b2 in itertools.combinations(term,2):
                    b3=[t for t in term if t not in (b1,b2)][0]
                    ps=even_paths(adj,N,b1,b2)
                    if ps is None: continue
                    if not ps: stat["no even path"]+=1; continue
                    avoid=[p for p in ps if b3 not in p]
                    # does the pair admit TWO internally disjoint even paths?
                    ok=False; okavoid=False
                    for i in range(len(ps)):
                        for j in range(i+1,len(ps)):
                            if set(ps[i][1:-1]) & set(ps[j][1:-1]): continue
                            ok=True
                            if b3 not in ps[i] and b3 not in ps[j]: okavoid=True; break
                        if okavoid: break
                    stat[(ok,okavoid)]+=1
                    if not ok:
                        detail.append((len(ps),len(avoid)))
    print("(two disjoint even paths exist?, and both avoiding b3?) -> count")
    for k,v in sorted(stat.items(),key=lambda x:str(x[0])): print("   ",k,"->",v)
    if detail:
        tot=[d[0] for d in detail]; av=[d[1] for d in detail]
        print("\nFAILING pairs: #even connecting paths  min=%d median=%d max=%d"%(min(tot),sorted(tot)[len(tot)//2],max(tot)))
        print("               of which AVOID b3       min=%d median=%d max=%d"%(min(av),sorted(av)[len(av)//2],max(av)))
        print("               fraction with ZERO avoiding b3: %d of %d"%(sum(1 for a in av if a==0),len(av)))

main()
