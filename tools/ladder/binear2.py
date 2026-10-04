"""CORRECTED Binary Ear Conjecture: allow a valuation SHIFT c.

If v_2(delta_j) = c+j-1 for j=1..k, then every delta_j is divisible by 2^c and
delta_j/2^c has valuation j-1, so by the Binary Ear Basis Lemma the subset sums
cover every multiple of 2^c modulo 2^{c+k}.  Hence, provided 2^c divides |C|,
exactly one subset sum L satisfies L = 0 mod 2^{c+k}; and if the window condition
   |C| + sum_{delta_j>0} delta_j < 2^{c+k+1}
holds, that L equals exactly 2^{c+k}.
(c=0 is the old statement; bipartite graphs need c>=1.)
"""
import itertools, sys

def v2(d):
    if d == 0: return None
    k = 0; d = abs(d)
    while d % 2 == 0: d //= 2; k += 1
    return k

def from_edges(n, E):
    adj=[[] for _ in range(n)]
    for a,b in E: adj[a].append(b); adj[b].append(a)
    return adj

def all_cycles(adj):
    n=len(adj); out=[]
    for s in range(n):
        st=[(s,[s],1<<s)]
        while st:
            u,path,mask=st.pop()
            for v in adj[u]:
                if v==s and len(path)>=3:
                    if path[1]<path[-1]: out.append(path[:])
                elif v>s and not (mask>>v)&1:
                    st.append((v,path+[v],mask|(1<<v)))
    return out

def ears_of(adj,cyc):
    n=len(adj); VC=set(cyc); L=len(cyc)
    ce=set()
    for i in range(L):
        a,b=cyc[i],cyc[(i+1)%L]; ce.add((min(a,b),max(a,b)))
    res=set()
    for u in cyc:
        for v in adj[u]:
            if v in VC and u<v and (u,v) not in ce: res.add((u,v,1,frozenset()))
    for u in cyc:
        for w in adj[u]:
            if w in VC: continue
            st=[(w,[w],{w})]
            while st:
                x,path,seen=st.pop()
                for y in adj[x]:
                    if y in VC:
                        if y!=u and u<y: res.add((u,y,len(path)+1,frozenset(path)))
                    elif y not in seen:
                        st.append((y,path+[y],seen|{y}))
    return list(res)

def arc_edges(cyc,u,v,which):
    L=len(cyc); pos={w:i for i,w in enumerate(cyc)}
    p,q=pos[u],pos[v]
    if which==0: idx=list(range(p,q)) if p<q else list(range(p,L))+list(range(0,q))
    else:        idx=list(range(q,p)) if q<p else list(range(q,L))+list(range(0,p))
    return frozenset((min(cyc[i],cyc[(i+1)%L]),max(cyc[i],cyc[(i+1)%L])) for i in idx)

def check(adj,name,kmax=6,cmax=4):
    best=None
    for cyc in all_cycles(adj):
        L=len(cyc)
        opts=[]
        for (u,v,ln,inter) in ears_of(adj,cyc):
            for which in (0,1):
                ae=arc_edges(cyc,u,v,which); d=ln-len(ae); val=v2(d)
                if val is None: continue
                opts.append((d,val,ae,inter))
        byval={}
        for o in opts: byval.setdefault(o[1],[]).append(o)
        for c in range(0,cmax+1):
            if L % (2**c) != 0: continue
            for k in range(2,kmax+1):
                need=[c+j for j in range(k)]
                if any(t not in byval for t in need): continue
                for combo in itertools.product(*[byval[t] for t in need]):
                    ok=True
                    for a,b in itertools.combinations(combo,2):
                        if a[2]&b[2] or a[3]&b[3]: ok=False;break
                    if not ok: continue
                    ps=sum(d for (d,_,_,_) in combo if d>0)
                    if L+ps < 2**(c+k+1):
                        tgt=2**(c+k)
                        if best is None or tgt>best[0]:
                            best=(tgt,c,k,L,[d for (d,_,_,_) in combo])
    if best:
        print("%-16s SATISFIED: target 2^%d=%d  (c=%d,k=%d, |C|=%d, deltas=%s)"
              % (name,(best[1]+best[2]),best[0],best[1],best[2],best[3],best[4]))
        return True
    print("%-16s *** NO BINARY EAR BASIS (any shift c) ***" % name)
    return False

def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return from_edges(n,sorted({tuple(sorted(e)) for e in E}))
def gp(n,k):
    E=[]
    for i in range(n): E+=[(i,(i+1)%n),(i,n+i),(n+i,n+((i+k)%n))]
    return from_edges(2*n,sorted({tuple(sorted(e)) for e in E}))

tests=[("K4",from_edges(4,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)])),
       ("K3,3",from_edges(6,[(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)])),
       ("Q3 (cube)",lcf(8,[3,-3])),
       ("PETERSEN",gp(5,2)),
       ("Heawood",lcf(14,[5,-5])),
       ("GP(7,2)",gp(7,2)),
       ("Mobius-Kantor",lcf(16,[5,-5])),
       ("Pappus",lcf(18,[5,7,-7,7,-7,-5])),
       ("Desargues",lcf(20,[5,-5,9,-9])),
       ("Dodecahedron",lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4]))]
print("CORRECTED BINARY EAR CONJECTURE -- exhaustive check\n")
bad=[]
for nm,g in tests:
    if not check(g,nm): bad.append(nm)
    sys.stdout.flush()
print()
print("failures:",bad)
