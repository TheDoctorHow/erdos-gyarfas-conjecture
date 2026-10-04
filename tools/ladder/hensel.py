"""HENSEL DESCENT test.

Move: from cycle C with ear P spanning arc A, get cycle C' of length |C|+delta,
delta = |P|-|A|.  If v_2(|C|-T) = v_2(delta) = t then v_2(|C'|-T) > t.

HYPOTHESIS (H_m): for every cycle C with 3 <= |C| < 2^{m+1} and |C| != 2^m there
is an ear P on C with v_2(delta) = v_2(|C| - 2^m) and 3 <= |C|+delta < 2^{m+1}.
H_m implies G has a cycle of length exactly 2^m.

NOTE: single-ear moves need NO disjointness between ears -- each move acts on a
fresh cycle.  That is the whole point of the revision.
"""
import itertools, random, sys

def v2(d):
    if d==0: return 10**9
    k=0; d=abs(d)
    while d%2==0: d//=2; k+=1
    return k
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
def ears_of(adj,cyc):
    VC=set(cyc); L=len(cyc); ce=set()
    for i in range(L):
        a,b=cyc[i],cyc[(i+1)%L]; ce.add((min(a,b),max(a,b)))
    res=set()
    for u in cyc:
        for v in adj[u]:
            if v in VC and u<v and (u,v) not in ce: res.add((u,v,1))
    for u in cyc:
        for w in adj[u]:
            if w in VC: continue
            st=[(w,[w],{w})]
            while st:
                x,path,seen=st.pop()
                for y in adj[x]:
                    if y in VC:
                        if y!=u and u<y: res.add((u,y,len(path)+1))
                    elif y not in seen: st.append((y,path+[y],seen|{y}))
    return list(res)
def deltas_on(adj,cyc):
    L=len(cyc); pos={w:i for i,w in enumerate(cyc)}
    out=[]
    for (u,v,ln) in ears_of(adj,cyc):
        p,q=pos[u],pos[v]
        a1=(q-p)%L; a2=L-a1
        out.append(ln-a1); out.append(ln-a2)
    return out

def test_Hm(adj,m):
    """check hypothesis H_m; return (holds, list of violating cycle lengths)"""
    T=2**m; hi=2**(m+1); viol=[]
    for cyc in all_cycles(adj):
        L=len(cyc)
        if L<3 or L>=hi or L==T: continue
        t=v2(L-T)
        ok=False
        for d in deltas_on(adj,cyc):
            if v2(d)==t and 3<=L+d<hi: ok=True; break
        if not ok: viol.append(L)
    return (len(viol)==0), sorted(set(viol))

def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return from_edges(n,sorted({tuple(sorted(e)) for e in E}))
def gp(n,k):
    E=[]
    for i in range(n): E+=[(i,(i+1)%n),(i,n+i),(n+i,n+((i+k)%n))]
    return from_edges(2*n,sorted({tuple(sorted(e)) for e in E}))

named=[("K4",from_edges(4,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)])),
       ("K3,3",from_edges(6,[(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)])),
       ("Q3 (cube)",lcf(8,[3,-3])),("PETERSEN",gp(5,2)),
       ("Heawood",lcf(14,[5,-5])),("GP(7,2)",gp(7,2)),
       ("Mobius-Kantor",lcf(16,[5,-5])),("Pappus",lcf(18,[5,7,-7,7,-7,-5])),
       ("Desargues",lcf(20,[5,-5,9,-9])),
       ("Dodecahedron",lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4]))]

print("HENSEL DESCENT -- hypothesis H_m per graph and target 2^m\n")
for nm,g in named:
    n=len(g); msgs=[]
    for m in range(2,7):
        if 2**m> n: break
        holds,viol=test_Hm(g,m)
        msgs.append("2^%d:%s%s"%(m,"OK" if holds else "FAIL",
                                 "" if holds else str(viol)))
    print("%-16s n=%-3d  %s"%(nm,n," | ".join(msgs)))
    sys.stdout.flush()
