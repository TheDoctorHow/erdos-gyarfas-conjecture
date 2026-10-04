import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools, sys
from collections import deque, Counter
base=_HERE+""
exec(open(base+"surv.py").read().split("def main")[0])
exec(open(base+"surv2.py").read().split("def main")[0].split("exec(")[0])

def build_B(nH,E,M):
    Ms=set(M); adj=[[] for _ in range(nH+len(M))]; term=[]; nxt=nH
    for (a,b) in E:
        if (a,b) in Ms:
            adj[a].append(nxt); adj[nxt].append(a)
            adj[b].append(nxt); adj[nxt].append(b); term.append(nxt); nxt+=1
        else:
            adj[a].append(b); adj[b].append(a)
    return adj,term

def even_paths(adj,n,s,t,limit=60000):
    out=[]; stack=[(s,1<<s,(s,))]
    while stack and len(out)<limit:
        u,mask,pth=stack.pop()
        for w in adj[u]:
            if w==t:
                if (len(pth))%2==0: out.append(pth+(t,))
            elif not (mask>>w & 1):
                stack.append((w,mask|1<<w,pth+(w,)))
    return out

def eset(p): return frozenset((min(p[i],p[i+1]),max(p[i],p[i+1])) for i in range(len(p)-1))

def comps(edges):
    """number of connected components of an edge set, and sizes"""
    adj={}
    for a,b in edges: adj.setdefault(a,[]).append(b); adj.setdefault(b,[]).append(a)
    seen=set(); sizes=[]
    for v in adj:
        if v in seen: continue
        q=[v]; seen.add(v); c=0
        while q:
            u=q.pop()
            for w in adj[u]:
                c+=1
                if w not in seen: seen.add(w); q.append(w)
        sizes.append(c//2)
    return sizes

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 4)
    stats=Counter(); segpar=Counter(); found=0
    for nH in (18,20):
        for trial in range(150):
            if found>=25: break
            g=randcubic(nH,rng)
            if not g: continue
            E,adj=g
            D=[dist(adj,nH,v) for v in range(nH)]
            def edist(e,f): return min(D[a][b] for a in e for b in f)
            cand=[t for t in itertools.combinations(E,3)
                  if all(edist(t[i],t[j])>=3 for i in range(3) for j in range(i+1,3))]
            if not cand: continue
            rng.shuffle(cand)
            for M in cand[:2]:
                if found>=25: break
                bad,term=build_B(nH,E,M)
                N=nH+len(M)
                for x,y in itertools.combinations(term,2):
                    ps=even_paths(bad,N,x,y)
                    if len(ps)<2: continue
                    byc={}
                    for p in ps: byc.setdefault((len(p)-1)%4,[]).append(p)
                    if len(byc)<2: continue
                    # minimal symmetric difference across classes
                    ks=sorted(byc)
                    bestd=None
                    for p in byc[ks[0]][:150]:
                        for qq in byc[ks[1]][:150]:
                            d=eset(p)^eset(qq)
                            if bestd is None or len(d)<len(bestd[0]): bestd=(d,p,qq)
                    d,p,qq=bestd
                    sizes=comps(d)
                    seg=len(eset(p)&d)
                    stats[(len(sizes),len(d)%4)]+=1
                    if len(sizes)==1:
                        segpar[(len(d)%4, seg%2)]+=1
                    found+=1
                    break
    print("minimal cross-class witnesses collected:",found)
    print("\n(#components of P xor P', |P xor P'| mod 4) -> count")
    for k,v in sorted(stats.items()): print("   ",k,"->",v)
    print("\nsingle-cycle witnesses:  (|Z| mod 4, segment length mod 2) -> count")
    for k,v in sorted(segpar.items()): print("   ",k,"->",v)

main()
