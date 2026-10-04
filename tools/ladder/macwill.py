import cmath

def gp(n,k):
    E=[]
    for i in range(n): E+=[(i,(i+1)%n),(i,n+i),(n+i,n+((i+k)%n))]
    return 2*n, sorted({tuple(sorted(e)) for e in E})
def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return n, sorted({tuple(sorted(e)) for e in E})

def W_direct(n,E,z):
    """sum over cycle space of z^{|x|}, by enumerating even subgraphs via a basis"""
    m=len(E); idx={e:i for i,e in enumerate(E)}
    adj={}
    for a,b in E: adj.setdefault(a,[]).append(b); adj.setdefault(b,[]).append(a)
    par={0:None}; pe={}; tree=set(); q=[0]
    while q:
        u=q.pop(0)
        for v in adj[u]:
            if v not in par:
                par[v]=u; e=tuple(sorted((u,v))); pe[v]=idx[e]; tree.add(idx[e]); q.append(v)
    def rootpath(v):
        s=set()
        while par[v] is not None: s.add(pe[v]); v=par[v]
        return s
    basis=[]
    for e in E:
        i=idx[e]
        if i in tree: continue
        a,b=e; cyc=(rootpath(a)^rootpath(b))|{i}
        msk=0
        for j in cyc: msk|=1<<j
        basis.append(msk)
    r=len(basis); tot=0
    for sub in range(1<<r):
        msk=0; s=sub; j=0
        while s:
            if s&1: msk^=basis[j]
            s>>=1; j+=1
        tot+=z**bin(msk).count("1")
    return tot,r

def W_via_cuts(n,E,z):
    """MacWilliams: W(z) = (1+z)^m / 2^{n-1} * sum over cuts of ((1-z)/(1+z))^{|cut|}"""
    m=len(E); ratio=(1-z)/(1+z); tot=0
    for S in range(1<<(n-1)):          # S ranges over subsets containing vertex 0 fixed out
        cut=0
        for a,b in E:
            ina=(S>>a)&1 if a<n-1 else 0
            inb=(S>>b)&1 if b<n-1 else 0
            if ina!=inb: cut+=1
        tot+=ratio**cut
    return ((1+z)**m)*tot/(2**(n-1))

z=cmath.exp(2j*cmath.pi/8)
print("ratio (1-z8)/(1+z8) =", (1-z)/(1+z), "   |.| =", abs((1-z)/(1+z)),
      "   tan(pi/8) =", cmath.tan(cmath.pi/8).real)
print()
for nm,(n,E) in [("Petersen",gp(5,2)),("Heawood",lcf(14,[5,-5])),("GP(7,2)",gp(7,2))]:
    for zz,lab in [(z,"zeta_8"),(1j,"i")]:
        a,r=W_direct(n,E,zz); b=W_via_cuts(n,E,zz)
        print("%-10s %-7s direct = %9.4f%+9.4fi   MacWilliams/cuts = %9.4f%+9.4fi   match=%s"
              % (nm,lab,a.real,a.imag,b.real,b.imag, abs(a-b)<1e-6*max(1,abs(a))))
    print()
