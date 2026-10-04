"""(H1*-short): does a SHALLOW cycle ever lack a negative shift?

Search many cubic 3-connected girth->=5 graphs.  For every shallow cycle record
whether min(shifts) < 0.  A shallow cycle with ALL shifts POSITIVE would be the
residual case; a shallow cycle with all shifts negative is already covered.
"""
import random, itertools, sys
from collections import Counter

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

random.seed(77)
graphs=0; shallow=0; nonegative=0; sizes=Counter()
for trial in range(200000):
    n=random.choice([14,16,18,20,22])
    perm=list(range(n)); random.shuffle(perm)
    pr=[(perm[2*i],perm[2*i+1]) for i in range(n//2)]
    bad=False
    for a,b in pr:
        s=min((b-a)%n,(a-b)%n)
        if s<4: bad=True; break
    if bad: continue
    if len(set(tuple(sorted(p)) for p in pr))!=n//2: continue
    adj=build(n,pr)
    if any(len(x)!=3 for x in adj): continue
    if girth(adj)<5: continue
    if not three_conn(adj): continue
    graphs+=1
    for cyc in all_cycles(adj):
        L=len(cyc)
        if ispow2(L): continue
        t=v2(L); ds=shifts(adj,cyc)
        if not ds: continue
        if max(v2(d) for d in ds) >= t: continue
        shallow+=1; sizes[(n,L,t)]+=1
        if min(ds) >= 0:
            nonegative+=1
            print("*** SHALLOW WITH NO NEGATIVE SHIFT *** n=%d |C|=%d t=%d shifts=%s"
                  %(n,L,t,sorted(set(ds))))
    if graphs>=120: break
print()
print("graphs tested:",graphs)
print("shallow cycles found:",shallow)
print("shallow cycles with NO negative shift:",nonegative)
print("shallow (n,|C|,t) profile:",dict(sizes))
