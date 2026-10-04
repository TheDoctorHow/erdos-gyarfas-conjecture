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

def shortest_par_path(adj,N,s,t,par,banned=frozenset()):
    """shortest s-t PATH of given parity, avoiding banned; BFS on parity-doubled graph
       (gives a walk; verified to be a path afterwards, else None)"""
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
    if len(set(path))!=len(path): return None      # walk, not a path
    return path

def is_2connected(adj,verts):
    verts=set(verts)
    if len(verts)<3: return False
    def conn(skip):
        st=[v for v in verts if v!=skip]
        if not st: return True
        seen={st[0]}; q=deque([st[0]])
        while q:
            u=q.popleft()
            for w in adj[u]:
                if w in verts and w!=skip and w not in seen: seen.add(w); q.append(w)
        return len(seen)==len(st)
    return all(conn(v) for v in list(verts)+[None])

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
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 88)
    st=Counter(); tot=0; goodpair=0; cfgok=0; cfgs=0
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
                cfgs+=1; anygood=False
                for b1,b2 in itertools.combinations(term,2):
                    tot+=1
                    alpha=shortest_par_path(adj,N,b1,b2,0)
                    if alpha is None: st["no even path at all"]+=1; continue
                    interior=set(alpha[1:-1])
                    rest=set(range(N))-interior
                    bip=is_bip(adj,rest); tc=is_2connected(adj,rest)
                    P=shortest_par_path(adj,N,b1,b2,0,banned=interior)
                    ok = P is not None and len(P)>1
                    st[("nonbip" if not bip else "BIPARTITE",
                        "2conn" if tc else "not2conn",
                        "even P exists" if ok else "NO even P")]+=1
                    if ok: goodpair+=1; anygood=True
                if anygood: cfgok+=1
    print("terminal pairs tested                    :",tot)
    print("  shortest even arc admits a DISJOINT even path (=> both arcs even):",goodpair)
    print("configurations tested                    :",cfgs)
    print("  with at least one such pair            :",cfgok)
    print("\n(bipartiteness of B_alpha, 2-connectivity of B_alpha, disjoint even path?) -> count")
    for k,v in sorted(st.items(),key=lambda x:-x[1]): print("   ",k,"->",v)

main()
