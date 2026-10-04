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

def escapes(adj,N,Zset,start,terms):
    seen={start}; q=deque([start])
    while q:
        u=q.popleft()
        for w in adj[u]:
            if w in Zset or w in seen: continue
            seen.add(w); q.append(w)
            if w in terms: return True
    return False

def cycles_through_pair(adj,N,b1,b2,cap=400000):
    """cycles containing both degree-2 vertices b1,b2; yields (cycleset, arc1len)"""
    x,y=adj[b1]; out=[]; cnt=0
    # a cycle through b1 is b1 + (x..y path); require it to contain b2
    stack=[(x,1<<x|1<<b1,(x,))]
    while stack:
        u,mask,pth=stack.pop()
        cnt+=1
        if cnt>cap: return out
        for w in adj[u]:
            if w==y:
                full=(b1,)+pth+(y,)
                if b2 in full: out.append(full)
            elif not (mask>>w & 1):
                stack.append((w,mask|1<<w,pth+(w,)))
    return out

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 33)
    tot=0; hit=0; miss=[]
    for nH in (18,20):
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
                tot+=1; found=False
                for b1,b2 in itertools.permutations(term,2):
                    if found: break
                    b3=[t for t in term if t not in (b1,b2)][0]
                    for cyc in cycles_through_pair(adj,N,b1,b2):
                        L=len(cyc)                    # cycle length in edges = #vertices
                        if L%2: continue
                        i1=cyc.index(b1); i2=cyc.index(b2)
                        a=(i2-i1)%L
                        if a%2 or (L-a)%2: continue   # need BOTH arcs even
                        Zset=set(cyc)
                        if b3 in Zset: continue
                        x,y=adj[b1]
                        if any(st not in Zset for st in (x,y)): pass
                        ok=False
                        for st in (x,y):
                            if st in Zset and escapes(adj,N,Zset,st,terms): ok=True; break
                        if ok:
                            found=True; break
                if found: hit+=1
                else: miss.append(nH)
    print("surviving configurations tested                      :",tot)
    print("  with a 2-terminal cycle, both arcs EVEN, escaping  :",hit)
    print("  without                                            :",len(miss))
    if miss: print("   misses by nH:",Counter(miss))

main()
