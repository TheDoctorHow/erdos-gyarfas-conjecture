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

def comps(adj,N,removed):
    seen=set(removed); out=[]
    for v in range(N):
        if v in seen: continue
        c={v}; q=deque([v]); seen.add(v)
        while q:
            u=q.popleft()
            for w in adj[u]:
                if w in seen: continue
                seen.add(w); c.add(w); q.append(w)
        out.append(c)
    return out

def bipartite_on(adj,N,verts):
    col={}
    for s in verts:
        if s in col: continue
        col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in adj[u]:
                if w not in verts: continue
                if w not in col: col[w]=1-col[u]; q.append(w)
                elif col[w]==col[u]: return False
    return True

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 55)
    tot=0; has3cut=0; hasfence=0; fence_termsplit=Counter()
    for nH in (18,20,22):
        for trial in range(120):
            g=randcubic(nH,rng)
            if not g: continue
            E,adj0=g
            D=[dist(adj0,nH,v) for v in range(nH)]
            def edist(e,f): return min(D[a][b] for a in e for b in f)
            cand=[t for t in itertools.combinations(E,3)
                  if all(edist(t[i],t[j])>=3 for i in range(3) for j in range(i+1,3))]
            if not cand: continue
            rng.shuffle(cand)
            for M in cand[:3]:
                adj,term=build_B(nH,E,M); N=nH+len(M); terms=set(term)
                tot+=1
                deg3=[v for v in range(N) if len(adj[v])==3]
                found3=False; foundfence=False
                for o,op in itertools.combinations(deg3,2):
                    cs=comps(adj,N,{o,op})
                    if len(cs)!=3: continue
                    found3=True
                    for K in cs:
                        rest=set(range(N))-K
                        if bipartite_on(adj,N,rest):
                            foundfence=True
                            # how are the three terminals distributed?
                            where=[]
                            for t in term:
                                for i,c in enumerate(cs):
                                    if t in c: where.append(i)
                            wi=[i for i,c in enumerate(cs) if c is K][0]
                            inodd=sum(1 for x in where if x==wi)
                            fence_termsplit[(len(set(where)),inodd)]+=1
                            break
                    if foundfence: break
                if found3: has3cut+=1
                if foundfence: hasfence+=1
    print("surviving configurations tested            :",tot)
    print("  with a 2-cut having THREE components     :",has3cut)
    print("  with a FENCE (3 comps + rest bipartite)  :",hasfence)
    print()
    print("for fences found: (#distinct components holding terminals, #terminals in odd part) -> count")
    for k,v in sorted(fence_termsplit.items()): print("   ",k,"->",v)

main()
