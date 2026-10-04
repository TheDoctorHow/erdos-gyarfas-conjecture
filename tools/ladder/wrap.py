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

def escapes(adj,Zset,start,terms):
    seen={start}; q=deque([start])
    while q:
        u=q.popleft()
        for w in adj[u]:
            if w in Zset or w in seen: continue
            seen.add(w); q.append(w)
            if w in terms: return True
    return False

def pair_ok(adj,b1,b2,terms,cap=300000):
    """is there a cycle through b1,b2 with BOTH arcs even, and a third terminal
       off it reachable?  returns (linkage_exists, linkage_and_escape)"""
    x,y=adj[b1]; cnt=0; link=False; full=False
    stack=[(x,1<<x|1<<b1,(x,))]
    while stack:
        u,mask,pth=stack.pop()
        cnt+=1
        if cnt>cap: break
        for w in adj[u]:
            if w==y:
                cyc=(b1,)+pth+(y,)
                if b2 not in cyc: continue
                L=len(cyc)
                i1=0; i2=cyc.index(b2)
                a=i2
                if a%2 or (L-a)%2: continue
                link=True
                Z=set(cyc)
                if any(t not in Z for t in terms):
                    for st in (x,y):
                        if st in Z and escapes(adj,Z,st,terms): full=True; break
                if full: return True,True
            elif not (mask>>w & 1):
                stack.append((w,mask|1<<w,pth+(w,)))
    return link,full

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 77)
    dist_link=Counter(); dist_full=Counter(); tot=0
    for nH in (18,20):
        for trial in range(100):
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
                adj,term=build_B(nH,E,M); terms=set(term)
                tot+=1; nl=0; nf=0
                for b1,b2 in itertools.combinations(term,2):
                    l1,f1=pair_ok(adj,b1,b2,terms)
                    l2,f2=pair_ok(adj,b2,b1,terms)
                    if l1 or l2: nl+=1
                    if f1 or f2: nf+=1
                dist_link[nl]+=1; dist_full[nf]+=1
    print("surviving configurations tested:",tot)
    print("\nnumber of the 3 terminal PAIRS admitting a both-arcs-even cycle:")
    for k in sorted(dist_link): print("   %d pairs -> %d configs"%(k,dist_link[k]))
    print("\nsame, additionally with the third-terminal ESCAPE:")
    for k in sorted(dist_full): print("   %d pairs -> %d configs"%(k,dist_full[k]))

main()
