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

def shortest_even(adj,N,s,t):
    seen={(s,0):None}; q=deque([(s,0)])
    while q:
        u,p=q.popleft()
        for w in adj[u]:
            st=(w,1-p)
            if st in seen: continue
            seen[st]=(u,p); q.append(st)
    if (t,0) not in seen: return None
    path=[]; cur=(t,0)
    while cur is not None: path.append(cur[0]); cur=seen[cur]
    path.reverse()
    return path if len(set(path))==len(path) else None

def two_disj_even(adj,N,s,t,cap=120000):
    ps=[]; stack=[(s,1<<s,(s,))]; c=0
    while stack:
        u,mask,p=stack.pop()
        c+=1
        if c>cap: break
        for w in adj[u]:
            if w==t:
                if len(p)%2==0: ps.append(p)
            elif not (mask>>w & 1): stack.append((w,mask|1<<w,p+(w,)))
    for i in range(len(ps)):
        for j in range(i+1,len(ps)):
            if not (set(ps[i][1:-1]) & set(ps[j][1:-1])): return True
    return False

def maxflow_to_set(adj,verts,src,T,k=3):
    """max number of internally disjoint paths from src to set T, capped at k
       (vertex-disjoint: split vertices)"""
    if src in T: return k
    idx={}; n=0
    for v in verts: idx[v]=n; n+=2
    cap={}
    def add(a,b,c):
        cap[(a,b)]=cap.get((a,b),0)+c; cap.setdefault((b,a),0)
    for v in verts:
        add(idx[v],idx[v]+1, k if (v==src or v in T) else 1)
    for v in verts:
        for w in adj[v]:
            if w in verts: add(idx[v]+1, idx[w], k)
    S=idx[src]+1; SINK=n
    for t in T: add(idx[t]+1, SINK, k)
    flow=0
    while flow<k:
        par={S:None}; q=deque([S]); found=False
        while q and not found:
            u=q.popleft()
            for (a,b),c in list(cap.items()):
                if a!=u or c<=0 or b in par: continue
                par[b]=u
                if b==SINK: found=True; break
                q.append(b)
        if not found: break
        cur=SINK
        while par[cur] is not None:
            p=par[cur]; cap[(p,cur)]-=1; cap[(cur,p)]=cap.get((cur,p),0)+1; cur=p
        flow+=1
    return flow

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 314)
    tested=0; failed=0; res=Counter()
    for nH in (18,20):
        for trial in range(160):
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
                for b1,b2 in itertools.permutations(term,2):
                    if two_disj_even(adj,N,b1,b2): continue      # only failing pairs
                    a=shortest_even(adj,N,b1,b2)
                    if a is None: continue
                    L=len(a)-1
                    inter=set(a[1:-1]); Bp=set(range(N))-inter
                    # T1 = b1 plus off-alpha neighbours of odd positions
                    T1={b1}
                    for i in range(1,L):
                        if i%2==1:
                            for u in adj[a[i]]:
                                if u not in a: T1.add(u)
                    tested+=1
                    bad=0
                    for w in Bp:
                        if w in T1: continue
                        if len(adj[w])!=3: continue
                        if maxflow_to_set(adj,Bp,w,T1,3)<3: bad+=1
                    res[bad]+=1
                    if bad: failed+=1
    print("failing pairs examined (so B' exists):",tested)
    print("  with SOME interior vertex lacking a fan to T1:",failed)
    print("\n#vertices lacking a fan -> count of instances")
    for k,v in sorted(res.items()): print("   ",k,"->",v)

main()
