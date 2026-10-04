import random, itertools, cmath

def girth(adj):
    n=len(adj); best=10**9
    for s in range(n):
        dist=[-1]*n; par=[-1]*n; dist[s]=0; q=[s]
        for u in q:
            for v in adj[u]:
                if dist[v]<0: dist[v]=dist[u]+1; par[v]=u; q.append(v)
                elif v!=par[u]: best=min(best,dist[u]+dist[v]+1)
    return best

def conn_after(adj,rem):
    n=len(adj); st0=next((v for v in range(n) if v not in rem),None)
    if st0 is None: return True
    seen={st0}; st=[st0]
    while st:
        u=st.pop()
        for v in adj[u]:
            if v not in rem and v not in seen: seen.add(v); st.append(v)
    return len(seen)==n-len(rem)

def three_conn(adj):
    n=len(adj)
    return n>=4 and all(conn_after(adj,set(p)) for p in itertools.combinations(range(n),2))

def bipartite(adj):
    n=len(adj); col=[-1]*n
    for s in range(n):
        if col[s]>=0: continue
        col[s]=0; q=[s]
        for u in q:
            for v in adj[u]:
                if col[v]<0: col[v]=1-col[u]; q.append(v)
                elif col[v]==col[u]: return False
    return True

def spectrum(adj):
    n=len(adj); c={}
    for s in range(n):
        stack=[(s,1,1<<s)]
        while stack:
            u,L,mask=stack.pop()
            for v in adj[u]:
                if v==s and L>=3: c[L]=c.get(L,0)+1
                elif v>s and not (mask>>v)&1: stack.append((v,L+1,mask|(1<<v)))
    return {L:k//2 for L,k in c.items()}

def from_edges(n,E):
    adj=[[] for _ in range(n)]
    for a,b in E: adj[a].append(b); adj[b].append(a)
    return adj

def gp(n,k):
    E=[]
    for i in range(n):
        E+= [(i,(i+1)%n),(i,n+i),(n+i,n+((i+k)%n))]
    return from_edges(2*n, sorted({tuple(sorted(e)) for e in E}))

def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return from_edges(n, sorted({tuple(sorted(e)) for e in E}))

def report(name,adj):
    spec=spectrum(adj); tot=sum(spec.values()); bip=bipartite(adj)
    # modulus 8 on |C|, or modulus 4 on |C|/2 when bipartite
    if bip:
        M=4; vals={L//2:c for L,c in spec.items()}
    else:
        M=8; vals=spec
    S=[]
    for k in range(1,M):
        s=sum(c*cmath.exp(2j*cmath.pi*k*L/M) for L,c in vals.items())
        S.append(abs(s)/tot)
    N0=sum(c for L,c in spec.items() if L%8==0)
    print("%-16s n=%-3d bip=%-5s M=%d  N0=%-6d  |S_k|/C: %s   max=%.3f  need<%.3f"
          % (name,len(adj),bip,M,N0," ".join("%.3f"%x for x in S),max(S),1.0/(M-1)))
    return max(S), 1.0/(M-1)

tests=[("Petersen",gp(5,2)),("Heawood",lcf(14,[5,-5])),("GP(7,2)",gp(7,2)),
       ("MobiusKantor",lcf(16,[5,-5])),("Pappus",lcf(18,[5,7,-7,7,-7,-5])),
       ("GP(9,2)",gp(9,2)),("Desargues",lcf(20,[5,-5,9,-9])),
       ("Dodecahedron",lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4])),
       ("GP(11,2)",gp(11,2)),("Nauru",lcf(24,[5,-9,7,-7,9,-5])),
       ("McGee",lcf(24,[12,7,-7]))]
print("=== character-sum profile (bipartite graphs reduced to |C|/2 mod 4) ===")
viol=0
for nm,a in tests:
    mx,thr=report(nm,a)
    if mx>=thr: viol+=1
random.seed(20260805)
for n in range(12,25,2):
    got=0;tries=0
    while got<3 and tries<8000:
        tries+=1
        perm=list(range(n)); random.shuffle(perm)
        pr=[(perm[2*i],perm[2*i+1]) for i in range(n//2)]
        if any((b-a)%n in (1,n-1) for a,b in pr): continue
        E=[(i,(i+1)%n) for i in range(n)]+[tuple(sorted(p)) for p in pr]
        if len(set(E))!=len(E): continue
        adj=from_edges(n,sorted(set(E)))
        if any(len(x)!=3 for x in adj): continue
        if girth(adj)<5 or not three_conn(adj): continue
        got+=1
        mx,thr=report("rand n=%d #%d"%(n,got),adj)
        if mx>=thr: viol+=1
print()
print("graphs violating max_k |S_k|/C < 1/(M-1):", viol)
