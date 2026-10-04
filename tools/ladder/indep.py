"""PREDICTION: a SHALLOW cycle must have G - V(C) an INDEPENDENT SET
(every bridge a single vertex), hence 3(n-|C|) <= |C|, i.e. |C| >= 3n/4.
Test on the four shallow Mobius-Kantor 12-cycles and on every shallow cycle found."""
import sys
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
def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return from_edges(n,sorted({tuple(sorted(e)) for e in E}))
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
        out += [ln-a1, ln-a2]
    return out

def test(adj,name):
    n=len(adj); found=0
    for cyc in all_cycles(adj):
        L=len(cyc)
        if ispow2(L): continue
        t=v2(L); ds=shifts(adj,cyc)
        if not ds: continue
        if max(v2(d) for d in ds) >= t: continue    # not shallow
        found+=1
        out=[v for v in range(n) if v not in set(cyc)]
        indep=all(v not in adj[u] for u in out for v in out)
        print("  %s shallow cycle |C|=%d t=%d  outside=%s  independent=%s  |C|>=3n/4? %s (3n/4=%.1f)"
              %(name,L,t,out,indep,L>=3*n/4,3*n/4))
        if found>=4: break
    if found==0: print("  %s: no shallow cycle"%name)

print("PREDICTION TEST: shallow => G-V(C) independent and |C| >= 3n/4\n")
for nm,g in [("Mobius-Kantor",lcf(16,[5,-5])),("Heawood",lcf(14,[5,-5])),
             ("Pappus",lcf(18,[5,7,-7,7,-7,-5])),("Desargues",lcf(20,[5,-5,9,-9])),
             ("Dodecahedron",lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4]))]:
    test(g,nm); sys.stdout.flush()
