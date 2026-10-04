"""IS (H1) TRUE?  (H1): a cycle whose complement contains an edge is never shallow.

Hand analysis suggests a bridge B = {b,b'} joined by an edge, with all FOUR of its
attachments in one position class mod 2^t, IS shallow:
   leaf ears  x1-b-x2   length 2, value 2 - 0 = 2 =/= 0
   cross ears x_i-b-b'-x_j length 3, value 3 - 0 = 3 =/= 0
So (H1) may be FALSE.  Search random cubic 3-connected girth->=5 graphs for a
SHALLOW cycle with a bridge of >= 2 vertices.
"""
import random, itertools, sys

def v2(d):
    if d==0: return 10**9
    k=0; d=abs(d)
    while d%2==0: d//=2; k+=1
    return k
def ispow2(x): return x>0 and (x&(x-1))==0
def girth(adj):
    n=len(adj); best=10**9
    for s in range(n):
        d=[-1]*n; p=[-1]*n; d[s]=0; q=[s]
        for u in q:
            for v in adj[u]:
                if d[v]<0: d[v]=d[u]+1; p[v]=u; q.append(v)
                elif v!=p[u]: best=min(best,d[u]+d[v]+1)
    return best
def conn_after(adj,rem):
    n=len(adj); s=next((v for v in range(n) if v not in rem),None)
    if s is None: return True
    seen={s}; st=[s]
    while st:
        u=st.pop()
        for v in adj[u]:
            if v not in rem and v not in seen: seen.add(v); st.append(v)
    return len(seen)==n-len(rem)
def three_conn(adj):
    return all(conn_after(adj,set(p)) for p in itertools.combinations(range(len(adj)),2))
def build(n,pairs):
    adj=[[] for _ in range(n)]
    for i in range(n):
        adj[i].append((i+1)%n); adj[(i+1)%n].append(i)
    for a,b in pairs: adj[a].append(b); adj[b].append(a)
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
def bridges_of(adj,cyc):
    """components of G - V(C); return list of (size, #attachments)"""
    VC=set(cyc); n=len(adj); seen=set(); res=[]
    for v in range(n):
        if v in VC or v in seen: continue
        comp={v}; st=[v]
        while st:
            u=st.pop()
            for w in adj[u]:
                if w not in VC and w not in comp: comp.add(w); st.append(w)
        seen|=comp
        att=sum(1 for u in comp for w in adj[u] if w in VC)
        res.append((len(comp),att))
    return res

random.seed(101)
print("SEARCH: shallow cycle with a bridge of >= 2 vertices  (would REFUTE (H1))\n")
found=0; graphs=0; shallowtot=0
for trial in range(300):
    n=random.choice([16,18,20,22,24])
    perm=list(range(n)); random.shuffle(perm)
    pr=[(perm[2*i],perm[2*i+1]) for i in range(n//2)]
    if any((b-a)%n in (1,n-1) for a,b in pr): continue
    if len(set(tuple(sorted(p)) for p in pr))!=n//2: continue
    adj=build(n,pr)
    if any(len(x)!=3 for x in adj): continue
    if girth(adj)<5 or not three_conn(adj): continue
    graphs+=1
    for cyc in all_cycles(adj):
        L=len(cyc)
        if ispow2(L): continue
        t=v2(L); ds=shifts(adj,cyc)
        if not ds: continue
        if max(v2(d) for d in ds) >= t: continue      # not shallow
        shallowtot+=1
        br=bridges_of(adj,cyc)
        big=[b for b in br if b[0]>=2]
        if big:
            found+=1
            print("*** SHALLOW with BIG BRIDGE ***  n=%d |C|=%d t=%d bridges(size,att)=%s"
                  %(n,L,t,br))
            print("    cycle:",cyc); print("    shifts:",sorted(set(ds)))
            if found>=5: break
    if found>=5: break
print()
print("graphs tested: %d   shallow cycles found: %d   with a bridge of >=2 vertices: %d"
      %(graphs,shallowtot,found))
