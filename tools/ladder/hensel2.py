"""HENSEL DESCENT, existential version.

H'_m : there EXISTS a cycle C_0 with 3 <= |C_0| < 2^{m+1} from which iterating
       'pick an ear with v_2(delta) = v_2(|C|-2^m), stay in the window'
       reaches a cycle of length exactly 2^m.
Each move strictly increases v_2(|C| - 2^m), so any descent has length <= m+1.
No disjointness is required: each move acts on a fresh cycle.
"""
import sys

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
def ears_with_paths(adj,cyc):
    """(u, v, [interior vertices in order from u to v])"""
    VC=set(cyc); L=len(cyc); ce=set()
    for i in range(L):
        a,b=cyc[i],cyc[(i+1)%L]; ce.add((min(a,b),max(a,b)))
    res=[]
    for u in cyc:
        for v in adj[u]:
            if v in VC and (min(u,v),max(u,v)) not in ce: res.append((u,v,[]))
    for u in cyc:
        for w in adj[u]:
            if w in VC: continue
            st=[(w,[w],{w})]
            while st:
                x,path,seen=st.pop()
                for y in adj[x]:
                    if y in VC:
                        if y!=u: res.append((u,y,path[:]))
                    elif y not in seen: st.append((y,path+[y],seen|{y}))
    return res

def reroute(cyc,u,v,inner,which):
    """replace one of the two u-v arcs of cyc by u-inner-v; return new cycle"""
    L=len(cyc); pos={w:i for i,w in enumerate(cyc)}
    p,q=pos[u],pos[v]
    # arc 'which'=0 : forward p->q ; keep the complementary part q->p
    if which==0:
        keep=[cyc[(q+i)%L] for i in range((p-q)%L+1)]     # v ... u
    else:
        keep=[cyc[(p+i)%L] for i in range((q-p)%L+1)]     # u ... v
        keep=keep[::-1]                                    # v ... u
    new=keep+inner[::-1] if False else keep+list(reversed(inner))
    # keep runs v...u ; then go u -> inner (reversed order gives u-side first?)
    new=keep+inner if keep[-1]==u else None
    if new is None: return None
    if len(set(new))!=len(new): return None
    return new

def descend(adj,m,verbose=False):
    T=2**m; hi=2**(m+1)
    for cyc in all_cycles(adj):
        if not (3<=len(cyc)<hi): continue
        cur=cyc[:]; path=[len(cur)]; ok=True
        for _ in range(m+2):
            L=len(cur)
            if L==T: return True,path
            t=v2(L-T); moved=False
            for (u,v,inner) in ears_with_paths(adj,cur):
                for which in (0,1):
                    nc=reroute(cur,u,v,inner,which)
                    if nc is None: continue
                    nl=len(nc)
                    d=nl-L
                    if v2(d)!=t: continue
                    if not (3<=nl<hi): continue
                    cur=nc; path.append(nl); moved=True; break
                if moved: break
            if not moved: ok=False; break
        if ok and len(cur)==T: return True,path
    return False,None

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
print("EXISTENTIAL HENSEL DESCENT: does SOME cycle descend to 2^m?\n")
for nm,g in named:
    res=[]
    for m in range(2,7):
        if 2**m>len(g): break
        ok,path=descend(g,m)
        res.append("2^%d:%s%s"%(m,"YES" if ok else "no", str(path) if ok else ""))
    print("%-16s %s"%(nm," | ".join(res)))
    sys.stdout.flush()
