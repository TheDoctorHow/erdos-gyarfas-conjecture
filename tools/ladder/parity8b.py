import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools

exec(open(_HERE+r"parity8.py".replace("parity8","parity8_lib")).read()) if False else None

def girth(adj):
    n = len(adj); best = 10**9
    for s in range(n):
        dist=[-1]*n; par=[-1]*n; dist[s]=0; q=[s]
        for u in q:
            for v in adj[u]:
                if dist[v]<0: dist[v]=dist[u]+1; par[v]=u; q.append(v)
                elif v!=par[u]: best=min(best,dist[u]+dist[v]+1)
    return best

def connected_after(adj, removed):
    n=len(adj); start=next((v for v in range(n) if v not in removed), None)
    if start is None: return True
    seen={start}; st=[start]
    while st:
        u=st.pop()
        for v in adj[u]:
            if v not in removed and v not in seen: seen.add(v); st.append(v)
    return len(seen)==n-len(removed)

def three_connected(adj):
    n=len(adj)
    if n<4: return False
    return all(connected_after(adj,set(p)) for p in itertools.combinations(range(n),2))

def cycles_by_edge(adj):
    """all cycles; return (N[8], per-edge N0 counts, per-vertex N0 counts)"""
    n=len(adj); N=[0]*8
    e0={}; v0={}
    for u in range(n):
        for w in adj[u]:
            if u<w: e0[(u,w)]=0
    for v in range(n): v0[v]=0
    for s in range(n):
        stack=[(s,[s],1<<s)]
        while stack:
            u,path,mask=stack.pop()
            for v in adj[u]:
                if v==s and len(path)>=3:
                    L=len(path); N[L%8]+=1
                    if L%8==0:
                        for i in range(len(path)):
                            a,b=path[i],path[(i+1)%len(path)]
                            e0[(min(a,b),max(a,b))]+=1
                        for a in path: v0[a]+=1
                elif v>s and not (mask>>v)&1:
                    stack.append((v,path+[v],mask|(1<<v)))
    N=[x//2 for x in N]
    e0={k:v//2 for k,v in e0.items()}
    v0={k:v//2 for k,v in v0.items()}
    return N,e0,v0

def from_edges(n,edges):
    adj=[[] for _ in range(n)]
    for a,b in edges: adj[a].append(b); adj[b].append(a)
    return adj

def gp(n,k):
    E=[]
    for i in range(n):
        E.append((i,(i+1)%n)); E.append((i,n+i)); E.append((n+i,n+((i+k)%n)))
    return from_edges(2*n, sorted({tuple(sorted(e)) for e in E}))

def lcf(n,shifts):
    E=[(i,(i+1)%n) for i in range(n)]
    for i in range(n): E.append((i,(i+shifts[i%len(shifts)])%n))
    return from_edges(n, sorted({tuple(sorted(e)) for e in E}))

def report(name,adj):
    g=girth(adj)
    N,e0,v0=cycles_by_edge(adj)
    tot=sum(N)
    mine=min(e0.values()); minv=min(v0.values())
    dev=max(abs(N[j]-tot/8) for j in range(8))/tot if tot else 0
    print("%-18s n=%-3d g=%-2d N0=%-6d par=%d  min_over_edges=%-5d min_over_verts=%-5d maxdev=%.3f"
          % (name,len(adj),g,N[0],N[0]%2,mine,minv,dev))
    return mine

print("=== edge-local and vertex-local 0-mod-8 cycle counts ===")
bad=0
tests=[("Petersen",gp(5,2)),("Heawood",lcf(14,[5,-5])),("GP(7,2)",gp(7,2)),
       ("MobiusKantor",lcf(16,[5,-5])),("Pappus",lcf(18,[5,7,-7,7,-7,-5])),
       ("GP(9,2)",gp(9,2)),("Desargues",lcf(20,[5,-5,9,-9])),
       ("Dodecahedron",lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4])),
       ("GP(11,2)",gp(11,2)),("Nauru",lcf(24,[5,-9,7,-7,9,-5])),
       ("McGee",lcf(24,[12,7,-7]))]
for nm,a in tests:
    if report(nm,a)==0: bad+=1

random.seed(20260805)
for n in range(10,25,2):
    got=0; tries=0
    while got<5 and tries<8000:
        tries+=1
        perm=list(range(n)); random.shuffle(perm)
        pairs=[(perm[2*i],perm[2*i+1]) for i in range(n//2)]
        if any((b-a)%n in (1,n-1) for a,b in pairs): continue
        E=[(i,(i+1)%n) for i in range(n)]+[tuple(sorted(p)) for p in pairs]
        if len(set(E))!=len(E): continue
        adj=from_edges(n,sorted(set(E)))
        if any(len(x)!=3 for x in adj): continue
        if girth(adj)<5 or not three_connected(adj): continue
        got+=1
        if report("rand n=%d #%d"%(n,got),adj)==0: bad+=1
print()
print("graphs with SOME edge on no 0-mod-8 cycle:", bad)
