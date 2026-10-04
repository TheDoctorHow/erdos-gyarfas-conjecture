import cmath, itertools, random

def from_edges(n,E):
    adj=[[] for _ in range(n)]
    for a,b in E: adj[a].append(b); adj[b].append(a)
    return adj
def gp(n,k):
    E=[]
    for i in range(n): E+=[(i,(i+1)%n),(i,n+i),(n+i,n+((i+k)%n))]
    return 2*n, sorted({tuple(sorted(e)) for e in E})
def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return n, sorted({tuple(sorted(e)) for e in E})

def analyse(name,n,E):
    m=len(E); idx={e:i for i,e in enumerate(E)}
    adj=from_edges(n,E)
    # spanning tree (BFS)
    par={0:None}; parent_edge={}; order=[0]; tree=set()
    while order:
        u=order.pop(0)
        for v in adj[u]:
            if v not in par:
                par[v]=u; e=tuple(sorted((u,v))); parent_edge[v]=idx[e]
                tree.add(idx[e]); order.append(v)
    # fundamental cycles for non-tree edges
    def path_to_root(v):
        s=set()
        while par[v] is not None:
            s.add(parent_edge[v]); v=par[v]
        return s
    basis=[]
    for e in E:
        i=idx[e]
        if i in tree: continue
        a,b=e
        pa,pb=path_to_root(a),path_to_root(b)
        cyc=(pa^pb)|{i}
        msk=0
        for j in cyc: msk|=1<<j
        basis.append(msk)
    r=len(basis); assert r==m-n+1
    # bicycle dimension: radical of Gram matrix G_ij = |B_i ^ B_j| mod 2
    G=[[bin(basis[i]&basis[j]).count("1")&1 for j in range(r)] for i in range(r)]
    # rank over F2
    rows=[sum(G[i][j]<<j for j in range(r)) for i in range(r)]
    rank=0; piv=[]
    for col in range(r):
        p=None
        for k in range(rank,r):
            if (rows[k]>>col)&1: p=k; break
        if p is None: continue
        rows[rank],rows[p]=rows[p],rows[rank]
        for k in range(r):
            if k!=rank and (rows[k]>>col)&1: rows[k]^=rows[rank]
        rank+=1
    d=r-rank   # bicycle dimension
    # enumerate the whole cycle space
    edges=E
    cnt={}          # (len mod 8, components) -> count
    lensum={}       # exact length -> count of elements
    bic_bad=False
    # recompute bicycle space vectors to test q on radical
    # radical basis from Gram
    # (rebuild: solve G v = 0)
    Gr=[sum(G[i][j]<<j for j in range(r)) for i in range(r)]
    rowsR=Gr[:]; rk=0; pivcols=[]
    for col in range(r):
        p=None
        for k in range(rk,r):
            if (rowsR[k]>>col)&1: p=k;break
        if p is None: continue
        rowsR[rk],rowsR[p]=rowsR[p],rowsR[rk]
        for k in range(r):
            if k!=rk and (rowsR[k]>>col)&1: rowsR[k]^=rowsR[rk]
        pivcols.append(col); rk+=1
    free=[c for c in range(r) if c not in pivcols]
    radbasis=[]
    for f in free:
        v=1<<f
        for i,c in enumerate(pivcols):
            if (rowsR[i]>>f)&1: v|=1<<c
        radbasis.append(v)
    # check q on radical: |x| mod 4 for all bicycles
    for sub in range(1<<len(radbasis)):
        msk=0
        for i in range(len(radbasis)):
            if (sub>>i)&1:
                sel=radbasis[i]
                for j in range(r):
                    if (sel>>j)&1: msk^=basis[j]
        if bin(msk).count("1")%4!=0: bic_bad=True
    # main enumeration
    for sub in range(1<<r):
        msk=0; s=sub; j=0
        while s:
            if s&1: msk^=basis[j]
            s>>=1; j+=1
        L=bin(msk).count("1")
        if L==0:
            cnt[(0,0)]=cnt.get((0,0),0)+1; continue
        # count components
        ad={}
        for i in range(m):
            if (msk>>i)&1:
                a,b=edges[i]; ad.setdefault(a,[]).append(b); ad.setdefault(b,[]).append(a)
        seen=set(); comp=0
        for v in ad:
            if v in seen: continue
            comp+=1; st=[v]; seen.add(v)
            while st:
                u=st.pop()
                for w in ad[u]:
                    if w not in seen: seen.add(w); st.append(w)
        cnt[(L%8,comp)]=cnt.get((L%8,comp),0)+1
        lensum[L]=lensum.get(L,0)+1
    tot=1<<r
    W=[sum(c*cmath.exp(2j*cmath.pi*k*L/8) for (L,comp),c in cnt.items()) for k in range(8)]
    maxcomp=max(comp for (L,comp) in cnt)
    S1=sum(c for (L,comp),c in cnt.items() if comp==1)
    N0=sum(c for (L,comp),c in cnt.items() if comp==1 and L==0)
    # exact mod-4 Gauss sum prediction
    pred = 0.0 if bic_bad else (2**d)*(2**((r-d)/2))
    print("%-14s n=%-3d r=%-3d |Z|=%-6d bicycle_dim=%d  q|rad!=0:%-5s"
          % (name,n,r,tot,d,bic_bad))
    print("      max components = %d  (bound floor(n/girth))   #cycles=%d  N0=%d" % (maxcomp,S1,N0))
    print("      |W(i)| computed = %.4f   predicted = %.4f" % (abs(W[2]),pred))
    print("      W(1)=%.0f  W(-1)=%.0f  |W(z8)|=%.3f |W(z8^2)|=%.3f |W(z8^3)|=%.3f"
          % (W[0].real,W[4].real,abs(W[1]),abs(W[2]),abs(W[3])))
    comps={}
    for (L,comp),c in cnt.items():
        if comp>0: comps[comp]=comps.get(comp,0)+c
    print("      elements by #components:", dict(sorted(comps.items())))
    print()

for nm,(n,E) in [("Petersen",gp(5,2)),("Heawood",lcf(14,[5,-5])),
                 ("GP(7,2)",gp(7,2)),("MobiusKantor",lcf(16,[5,-5])),
                 ("Pappus",lcf(18,[5,7,-7,7,-7,-5])),("GP(9,2)",gp(9,2)),
                 ("Dodecahedron",lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4])),
                 ("Desargues",lcf(20,[5,-5,9,-9]))]:
    analyse(nm,n,E)
