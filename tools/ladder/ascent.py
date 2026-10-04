"""ASCENT STEP census.

STEP: for a cycle C with |C| not a power of two and t = v_2(|C|), find an ear P
and arc A with v_2(|P|-|A|) = t.  (Sufficient: 2^t | |P| and 2^t | |A|.)

Census: over all cycles of each test graph, how often does the step FAIL, and what
do the failures look like?
"""
import sys, itertools
from collections import Counter

def v2(d):
    if d==0: return 10**9
    k=0; d=abs(d)
    while d%2==0: d//=2; k+=1
    return k
def ispow2(x): return x>0 and (x&(x-1))==0
def from_edges(n,E):
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
                elif v>s and not (mask>>v)&1: st.append((v,path+[v],mask|(1<<v)))
    return out
def shifts(adj,cyc):
    """multiset of available delta = |P|-|A| over all ears and both arcs"""
    VC=set(cyc); L=len(cyc); pos={w:i for i,w in enumerate(cyc)}
    ce=set()
    for i in range(L):
        a,b=cyc[i],cyc[(i+1)%L]; ce.add((min(a,b),max(a,b)))
    ears=[]
    for u in cyc:
        for v in adj[u]:
            if v in VC and u<v and (min(u,v),max(u,v)) not in ce: ears.append((u,v,1))
    for u in cyc:
        for w in adj[u]:
            if w in VC: continue
            st=[(w,[w],{w})]
            while st:
                x,path,seen=st.pop()
                for y in adj[x]:
                    if y in VC:
                        if y!=u and u<y: ears.append((u,y,len(path)+1))
                    elif y not in seen: st.append((y,path+[y],seen|{y}))
    out=[]
    for (u,v,ln) in ears:
        p,q=pos[u],pos[v]; a1=(q-p)%L; a2=L-a1
        out.append(ln-a1); out.append(ln-a2)
    return out

def census(adj,name):
    tot=0; fail=0; failinfo=Counter(); failex=[]
    for cyc in all_cycles(adj):
        L=len(cyc)
        if ispow2(L): continue
        tot+=1; t=v2(L)
        ds=shifts(adj,cyc)
        if any(v2(d)==t for d in ds): continue
        fail+=1
        vs=tuple(sorted(set(v2(d) for d in ds if d!=0)))
        failinfo[(L,t,vs)]+=1
        if len(failex)<3: failex.append((L,t,sorted(set(ds))))
    print("%-16s cycles(non-2^k)=%-5d  step FAILS on %d" % (name,tot,fail))
    for k,v in sorted(failinfo.items()):
        print("      |C|=%-3d t=%-2d available v_2 = %-12s  count=%d"%(k[0],k[1],str(k[2]),v))
    for (L,t,ds) in failex[:2]:
        print("      example |C|=%d t=%d shifts=%s"%(L,t,ds))
    return tot,fail

def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return from_edges(n,sorted({tuple(sorted(e)) for e in E}))
def gp(n,k):
    E=[]
    for i in range(n): E+=[(i,(i+1)%n),(i,n+i),(n+i,n+((i+k)%n))]
    return from_edges(2*n,sorted({tuple(sorted(e)) for e in E}))

print("ASCENT STEP census  (t = v_2(|C|); need an ear with v_2(delta) = t)\n")
for nm,g in [("K4",from_edges(4,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)])),
             ("K3,3",from_edges(6,[(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)])),
             ("Q3",lcf(8,[3,-3])),("Petersen",gp(5,2)),("Heawood",lcf(14,[5,-5])),
             ("GP(7,2)",gp(7,2)),("Mobius-Kantor",lcf(16,[5,-5])),
             ("Pappus",lcf(18,[5,7,-7,7,-7,-5])),("Desargues",lcf(20,[5,-5,9,-9])),
             ("Dodecahedron",lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4]))]:
    census(g,nm); sys.stdout.flush()
