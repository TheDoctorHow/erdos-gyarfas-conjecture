"""Binary Ear Conjecture, FIXED checker.
Fixes over binear2.py:
  (1) each EAR may be used at most once (previously the same chord could be taken
      twice, once with each of its two arcs, producing a degenerate 2-cycle);
  (2) every subset sum must give a length >= 3, i.e. |C| + sum of negative deltas >= 3
      (otherwise some of the 2^k 'cycles' do not exist and the bijection argument
      does not apply).
"""
import itertools, random, sys

def v2(d):
    if d==0: return None
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
                    elif y not in seen: st.append((y,path+[y],seen|{y}))
    return list(res)
def arc_edges(cyc,u,v,which):
    L=len(cyc); pos={w:i for i,w in enumerate(cyc)}
    p,q=pos[u],pos[v]
    if which==0: idx=list(range(p,q)) if p<q else list(range(p,L))+list(range(0,q))
    else:        idx=list(range(q,p)) if q<p else list(range(q,L))+list(range(0,p))
    return frozenset((min(cyc[i],cyc[(i+1)%L]),max(cyc[i],cyc[(i+1)%L])) for i in idx)

def check(adj,name,kmax=4,cmax=3,quiet=False):
    best=None
    for cyc in all_cycles(adj):
        L=len(cyc); opts=[]
        for eid,(u,v,ln,inter) in enumerate(ears_of(adj,cyc)):
            for which in (0,1):
                ae=arc_edges(cyc,u,v,which); d=ln-len(ae); val=v2(d)
                if val is None: continue
                opts.append((d,val,ae,inter,eid))
        byval={}
        for o in opts: byval.setdefault(o[1],[]).append(o)
        for c in range(0,cmax+1):
            if L % (2**c): continue
            for k in range(2,kmax+1):
                need=[c+j for j in range(k)]
                if any(t not in byval for t in need): continue
                for combo in itertools.product(*[byval[t] for t in need]):
                    if len({o[4] for o in combo})!=k: continue      # FIX (1)
                    ok=True
                    for a,b in itertools.combinations(combo,2):
                        if a[2]&b[2] or a[3]&b[3]: ok=False;break
                    if not ok: continue
                    ps=sum(d for (d,_,_,_,_) in combo if d>0)
                    ng=sum(d for (d,_,_,_,_) in combo if d<0)
                    if L+ng < 3: continue                            # FIX (2)
                    if L+ps < 2**(c+k+1):
                        tgt=2**(c+k)
                        if best is None or tgt>best[0]:
                            best=(tgt,c,k,L,[d for (d,_,_,_,_) in combo])
    if best and not quiet:
        print("%-16s SATISFIED: target %d  (c=%d,k=%d,|C|=%d,deltas=%s)"
              %(name,best[0],best[1],best[2],best[3],best[4]))
    elif not quiet:
        print("%-16s *** NO BINARY EAR BASIS ***"%name)
    return best is not None

def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return from_edges(n,sorted({tuple(sorted(e)) for e in E}))
def gp(n,k):
    E=[]
    for i in range(n): E+=[(i,(i+1)%n),(i,n+i),(n+i,n+((i+k)%n))]
    return from_edges(2*n,sorted({tuple(sorted(e)) for e in E}))

print("FIXED checker -- named graphs\n")
named=[("K4",from_edges(4,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)])),
       ("K3,3",from_edges(6,[(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)])),
       ("Q3 (cube)",lcf(8,[3,-3])),("PETERSEN",gp(5,2)),
       ("Heawood",lcf(14,[5,-5])),("GP(7,2)",gp(7,2)),
       ("Mobius-Kantor",lcf(16,[5,-5])),("Pappus",lcf(18,[5,7,-7,7,-7,-5]))]
bad=[]
for nm,g in named:
    if not check(g,nm): bad.append(nm)
    sys.stdout.flush()
print("\nnamed failures:",bad)

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

print("\nFIXED checker -- random 3-connected cubic girth-5 graphs\n")
random.seed(3); fails=[]
for n in range(12,21,2):
    got=0; tries=0
    while got<3 and tries<40000:
        tries+=1
        perm=list(range(n)); random.shuffle(perm)
        pr=[(perm[2*i],perm[2*i+1]) for i in range(n//2)]
        if any((b-a)%n in (1,n-1) for a,b in pr): continue
        if len(set(tuple(sorted(p)) for p in pr))!=n//2: continue
        adj=build(n,pr)
        if any(len(x)!=3 for x in adj): continue
        if girth(adj)<5 or not three_conn(adj): continue
        got+=1
        if not check(adj,"rand n=%d #%d"%(n,got)): fails.append((n,got))
        sys.stdout.flush()
print("\nrandom failures:",fails)
