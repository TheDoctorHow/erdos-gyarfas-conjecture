import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools, sys
from collections import deque, Counter
base=_HERE+""
exec(open(base+"surv.py").read().split("def main")[0])
exec(open(base+"place.py").read().split("def main")[0].split("exec(")[0])

def build_B(nH,E,M):
    Ms=set(M); adj=[[] for _ in range(nH+len(M))]; term=[]; nxt=nH
    for (a,b) in E:
        if (a,b) in Ms:
            adj[a].append(nxt); adj[nxt].append(a)
            adj[b].append(nxt); adj[nxt].append(b); term.append(nxt); nxt+=1
        else:
            adj[a].append(b); adj[b].append(a)
    return adj,term

def shortest_par_path(adj,N,s,t,par,banned=frozenset()):
    start=(s,0); seen={start:None}; q=deque([start])
    while q:
        (u,p)=q.popleft()
        for w in adj[u]:
            if w in banned: continue
            st=(w,1-p)
            if st in seen: continue
            seen[st]=(u,p); q.append(st)
    goal=(t,par)
    if goal not in seen: return None
    path=[]; cur=goal
    while cur is not None: path.append(cur[0]); cur=seen[cur]
    path.reverse()
    if len(set(path))!=len(path): return None
    return path

def is_bip(adj,verts):
    verts=set(verts); col={}
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
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 91)
    joint=Counter(); lens={True:[],False:[]}
    allthree=0; cfgs=0
    for nH in (18,20):
        for trial in range(140):
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
                cfgs+=1; fails=0
                for b1,b2 in itertools.combinations(term,2):
                    b3=[t for t in term if t not in (b1,b2)][0]
                    a=shortest_par_path(adj,N,b1,b2,0)
                    if a is None: continue
                    inter=set(a[1:-1])
                    isT = is_bip(adj,set(range(N))-inter)      # interior is a transversal
                    hasb3 = b3 in inter
                    joint[(isT,hasb3)]+=1
                    lens[isT].append(len(a)-1)
                    if isT: fails+=1
                if fails==3: allthree+=1
    print("configurations tested            :",cfgs)
    print("  with ALL THREE pairs failing   :",allthree)
    print()
    print("(interior is an odd-cycle transversal?, third terminal inside the interior?) -> count")
    for k,v in sorted(joint.items()): print("   ",k,"->",v)
    for k in (True,False):
        L=sorted(lens[k])
        if L: print("  |alpha| when transversal=%s : min=%d median=%d max=%d  (n=%d)"%(
            k,L[0],L[len(L)//2],L[-1],len(L)))

main()
