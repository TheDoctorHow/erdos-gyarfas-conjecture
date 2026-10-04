"""Diagnose the Mobius-Kantor obstruction, then test the WINDOWLESS criterion.

CLAIM (to be checked): the window [3,2^{m+1}) is NOT needed for correctness.
Termination only needs |C| <= n: if v_2(|C|-T) = t and |C| != T then 2^t <= |C|+T,
so t <= log_2(n+T); the valuation cannot increase past that without |C| = T.
So dropping the window only ENLARGES the admissible move set.
"""
import sys
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
def canon(c):
    n=len(c); best=None
    for r in (list(c), list(c)[::-1]):
        for i in range(n):
            tt=tuple(r[i:]+r[:i])
            if best is None or tt<best: best=tt
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

def reach(adj,m,use_window):
    T=2**m; hi=2**(m+1); n=len(adj)
    cycles=[c for c in all_cycles(adj) if (3<=len(c)<hi) or not use_window]
    idx={canon(c):i for i,c in enumerate(cycles)}
    rev=[[] for _ in cycles]
    for i,c in enumerate(cycles):
        t=v2(len(c)-T)
        for nc in neighbours(adj,c):
            if use_window and not (3<=len(nc)<hi): continue
            k=idx.get(canon(nc))
            if k is None: continue
            if v2(len(nc)-T)>=t: rev[k].append(i)
    tg=[i for i,c in enumerate(cycles) if len(c)==T]
    if not tg: return None
    seen=set(tg); dq=deque(tg)
    while dq:
        u=dq.popleft()
        for w in rev[u]:
            if w not in seen: seen.add(w); dq.append(w)
    bad=[cycles[i] for i in range(len(cycles)) if i not in seen]
    return len(cycles),len(seen),bad

def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return from_edges(n,sorted({tuple(sorted(e)) for e in E}))
MK=lcf(16,[5,-5])

print("=== DIAGNOSIS: Mobius-Kantor, T=8, WITH window [3,16) ===")
r=reach(MK,3,True); tot,rc,bad=r
print("cycles %d, reach %d, stuck %d  lengths %s"%(tot,rc,len(bad),sorted(len(c) for c in bad)))
for c in bad[:1]:
    print("  example stuck cycle:",c)
    t=v2(len(c)-8)
    print("  |C|=%d, t=v_2(|C|-8)=%d"%(len(c),t))
    ds={}
    for nc in neighbours(MK,c):
        d=len(nc)-len(c); ds.setdefault(d,0); ds[d]+=1
    print("  available shifts delta -> new length (count):")
    for d in sorted(ds):
        print("     delta=%-4d new len=%-3d  v_2=%-3s  in window=%s  count=%d"
              %(d,len(c)+d,("inf" if d==0 else v2(d)), 3<=len(c)+d<16, ds[d]))

print()
print("=== WINDOWLESS criterion (only 3 <= |C| <= n) ===")
for nm,g,ms in [("Mobius-Kantor",MK,[3,4]),
                ("Heawood",lcf(14,[5,-5]),[3]),
                ("Pappus",lcf(18,[5,7,-7,7,-7,-5]),[3,4]),
                ("Desargues",lcf(20,[5,-5,9,-9]),[3,4]),
                ("Dodecahedron",lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4]),[3,4])]:
    outs=[]
    for m in ms:
        r=reach(g,m,False)
        if r is None: outs.append("2^%d:no target"%m); continue
        tot,rc,bad=r
        outs.append("2^%d:%s(%d/%d)%s"%(m,"ALL" if not bad else "part",rc,tot,
                    "" if not bad else " stuck lengths "+str(sorted(set(len(c) for c in bad)))))
    print("%-16s %s"%(nm," | ".join(outs)))
    sys.stdout.flush()
