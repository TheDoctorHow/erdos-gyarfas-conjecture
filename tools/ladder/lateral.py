"""EXISTENTIAL STEP with LATERAL moves.

For target T = 2^m and window W = [3, 2^{m+1}), and a cycle C with t = v_2(|C|-T):
   PROGRESS ear : v_2(delta) = t      -> new level > t
   LATERAL  ear : v_2(delta) > t      -> new level = t   (free exploration)
   REGRESSIVE   : v_2(delta) < t      -> forbidden
A move is ADMISSIBLE if it is progress or lateral, i.e. iff t(C') >= t(C), and the
new length stays in W.

TEST: is the reachability UNIVERSAL, i.e. can EVERY cycle in W reach a cycle of
length exactly T using only admissible moves?  (Progress-only reachability was
already 26/26 from some start; the question is whether lateral moves make it work
from EVERY start, which would let us drop the 'choose C_0' freedom entirely.)
"""
import sys, random, itertools
from collections import deque

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
def canon(cyc):
    n=len(cyc); best=None
    for r in (cyc, cyc[::-1]):
        for i in range(n):
            c=tuple(r[i:]+r[:i])
            if best is None or c<best: best=c
    return best
def neighbours(adj,cyc):
    VC=set(cyc); L=len(cyc); pos={w:i for i,w in enumerate(cyc)}
    ce=set()
    for i in range(L):
        a,b=cyc[i],cyc[(i+1)%L]; ce.add((min(a,b),max(a,b)))
    ears=[]
    for u in cyc:
        for v in adj[u]:
            if v in VC and u<v and (min(u,v),max(u,v)) not in ce: ears.append((u,v,[]))
    for u in cyc:
        for w in adj[u]:
            if w in VC: continue
            st=[(w,[w],{w})]
            while st:
                x,path,seen=st.pop()
                for y in adj[x]:
                    if y in VC:
                        if y!=u: ears.append((u,y,path[:]))
                    elif y not in seen: st.append((y,path+[y],seen|{y}))
    out=[]
    for (u,v,inner) in ears:
        p,q=pos[u],pos[v]
        for which in (0,1):
            if which==0: keep=[cyc[(q+i)%L] for i in range((p-q)%L+1)]
            else:        keep=[cyc[(p+i)%L] for i in range((q-p)%L+1)][::-1]
            if keep[0]!=v or keep[-1]!=u: continue
            new=keep+inner
            if len(set(new))!=len(new) or len(new)<3: continue
            out.append(new)
    return out

def universal_reach(adj,m):
    T=2**m; hi=2**(m+1)
    cycles=[c for c in all_cycles(adj) if 3<=len(c)<hi]
    if not cycles: return None
    idx={canon(c):i for i,c in enumerate(cycles)}
    # build reverse edges: C -> C' admissible iff t(C') >= t(C)
    rev=[[] for _ in cycles]
    for i,c in enumerate(cycles):
        t=v2(len(c)-T)
        for nc in neighbours(adj,c):
            if not (3<=len(nc)<hi): continue
            k=idx.get(canon(nc))
            if k is None: continue
            if v2(len(nc)-T)>=t: rev[k].append(i)
    targets=[i for i,c in enumerate(cycles) if len(c)==T]
    if not targets: return (False,len(cycles),0,None)
    seen=set(targets); dq=deque(targets)
    while dq:
        u=dq.popleft()
        for w in rev[u]:
            if w not in seen: seen.add(w); dq.append(w)
    bad=[len(cycles[i]) for i in range(len(cycles)) if i not in seen]
    return (len(bad)==0, len(cycles), len(seen), sorted(set(bad))[:8])

def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return from_edges(n,sorted({tuple(sorted(e)) for e in E}))
def gp(n,k):
    E=[]
    for i in range(n): E+=[(i,(i+1)%n),(i,n+i),(n+i,n+((i+k)%n))]
    return from_edges(2*n,sorted({tuple(sorted(e)) for e in E}))

named=[("K4",from_edges(4,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)])),
       ("K3,3",from_edges(6,[(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)])),
       ("Q3",lcf(8,[3,-3])),("Petersen",gp(5,2)),("Heawood",lcf(14,[5,-5])),
       ("GP(7,2)",gp(7,2)),("Mobius-Kantor",lcf(16,[5,-5])),
       ("Pappus",lcf(18,[5,7,-7,7,-7,-5])),("Desargues",lcf(20,[5,-5,9,-9])),
       ("Dodecahedron",lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4]))]
print("UNIVERSAL reachability of 2^m under admissible (progress OR lateral) moves\n")
for nm,g in named:
    lines=[]
    for m in range(2,6):
        if 2**m>len(g): break
        r=universal_reach(g,m)
        if r is None: continue
        ok,tot,reach,bad=r
        lines.append("2^%d:%s(%d/%d)%s"%(m,"ALL" if ok else "part",reach,tot,
                                         "" if ok else " missing from lengths "+str(bad)))
    print("%-16s %s"%(nm," | ".join(lines)))
    sys.stdout.flush()
