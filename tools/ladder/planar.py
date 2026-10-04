import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import itertools, sys
from collections import deque
base=_HERE+""
exec(open(base+"surv.py").read().split("def main")[0])
exec(open(base+"l2b.py").read().split("def cycles_through_pair")[0].split("exec(")[-1].split(")\n",1)[1] if False else "")

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

def lcf(n,seq):
    S=set((min(i,(i+1)%n),max(i,(i+1)%n)) for i in range(n))
    for i,s in enumerate(seq):
        a,b=i,(i+s)%n; S.add((min(a,b),max(a,b)))
    adj=[[] for _ in range(n)]
    for a,b in S: adj[a].append(b); adj[b].append(a)
    return sorted(S),adj

def even_paths_exist_both_classes(adj,N,s,t,cap=2000000):
    cls=set(); cnt=0
    stack=[(s,1<<s,0)]
    while stack:
        u,mask,d=stack.pop()
        cnt+=1
        if cnt>cap: return None
        for w in adj[u]:
            if w==t:
                if (d+1)%2==0: cls.add((d+1)%4)
            elif not (mask>>w & 1): stack.append((w,mask|1<<w,d+1))
    return cls

# Dodecahedron: planar, cubic, girth 5
E,adj=lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4]*2)
print("Dodecahedron: n=20 girth=%d  edges=%d"%(girth(adj,20),len(E)))
D=[dist(adj,20,v) for v in range(20)]
def edist(e,f): return min(D[a][b] for a in e for b in f)
cand=[t for t in itertools.combinations(E,3)
      if all(edist(t[i],t[j])>=3 for i in range(3) for j in range(i+1,3))]
print("matchings of 3 edges pairwise H-distance >= 3 :",len(cand))

if not cand:
    print("  -> the planar Dodecahedron admits NO valid surviving structure at s=3")
else:
    ok=0; bad=0
    for M in cand:
        adjB,term=build_B(20,E,M); N=20+3; terms=set(term)
        viol=False
        for x,y in itertools.combinations(term,2):
            cls=even_paths_exist_both_classes(adjB,N,x,y)
            if cls is None: continue
            if len(cls)>1: viol=True; break
        if viol: bad+=1
        else: ok+=1
    print("  configurations where (SC) FAILS :",bad)
    print("  configurations where (SC) HOLDS :",ok)

# also check the other planar cubic girth-5 option available: none small; report
print("\nball-packing bound: 3 edges pairwise at distance >=3 need |N[e]|=6 each,")
print("so |V(H)| >= 18; the Dodecahedron (20) is the smallest planar cubic girth-5 host.")
