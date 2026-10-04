"""Hensel Descent Conjecture -- proper search (backtracking, not greedy).

MOVE SET (proved in 172/173): from cycle C, one ear-reroute reaches exactly the
cycles C' with C ^ C' a single path -- i.e. |C'| = |P|+|A| or |P|+|A'| for an ear
P spanning arcs A, A'.  A move is ADMISSIBLE for target T if
v_2(|C'| - T) > v_2(|C| - T) and 3 <= |C'| < 2T.
Since v_2 strictly increases and is <= m while |C| != T, depth <= m+1.
"""
import sys, random, itertools

def v2(d):
    if d == 0: return 10**9
    k = 0; d = abs(d)
    while d % 2 == 0: d //= 2; k += 1
    return k
def from_edges(n, E):
    adj = [[] for _ in range(n)]
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

def neighbours(adj, cyc):
    """all cycles reachable by one ear-reroute; returns list of cycles"""
    VC=set(cyc); L=len(cyc); pos={w:i for i,w in enumerate(cyc)}
    ce=set()
    for i in range(L):
        a,b=cyc[i],cyc[(i+1)%L]; ce.add((min(a,b),max(a,b)))
    ears=[]
    for u in cyc:                                    # chords
        for v in adj[u]:
            if v in VC and u<v and (min(u,v),max(u,v)) not in ce:
                ears.append((u,v,[]))
    for u in cyc:                                    # ears through outside
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
            # keep the arc from v round to u (which=0) or from u round to v then reverse
            if which==0: keep=[cyc[(q+i)%L] for i in range((p-q)%L+1)]   # v..u
            else:        keep=[cyc[(p+i)%L] for i in range((q-p)%L+1)][::-1]  # v..u
            if keep[0]!=v or keep[-1]!=u: continue
            new=keep+inner
            if len(set(new))!=len(new): continue
            if len(new)<3: continue
            out.append(new)
    return out

def descend(adj, m, maxnodes=200000):
    T=2**m; hi=2**(m+1); nodes=0
    def dfs(cyc, depth):
        nonlocal nodes
        nodes+=1
        if nodes>maxnodes: return None
        L=len(cyc)
        if L==T: return [L]
        if depth<=0: return None
        t=v2(L-T)
        for nc in neighbours(adj,cyc):
            nl=len(nc)
            if not (3<=nl<hi): continue
            if v2(nl-T)<=t: continue
            r=dfs(nc, depth-1)
            if r is not None: return [L]+r
        return None
    for cyc in all_cycles(adj):
        if not (3<=len(cyc)<hi): continue
        r=dfs(cyc, m+1)
        if r is not None: return True, r
    return False, None

def lcf(n,sh):
    E=[(i,(i+1)%n) for i in range(n)]+[(i,(i+sh[i%len(sh)])%n) for i in range(n)]
    return from_edges(n,sorted({tuple(sorted(e)) for e in E}))
def gp(n,k):
    E=[]
    for i in range(n): E+=[(i,(i+1)%n),(i,n+i),(n+i,n+((i+k)%n))]
    return from_edges(2*n,sorted({tuple(sorted(e)) for e in E}))
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

named=[("K4",from_edges(4,[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)])),
       ("K3,3",from_edges(6,[(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)])),
       ("Q3",lcf(8,[3,-3])),("Petersen",gp(5,2)),("Heawood",lcf(14,[5,-5])),
       ("GP(7,2)",gp(7,2)),("Mobius-Kantor",lcf(16,[5,-5])),
       ("Pappus",lcf(18,[5,7,-7,7,-7,-5])),("Desargues",lcf(20,[5,-5,9,-9])),
       ("Dodecahedron",lcf(20,[10,7,4,-4,-7,10,-4,7,-7,4])),("McGee",lcf(24,[12,7,-7]))]
print("HENSEL DESCENT with backtracking -- named graphs\n")
for nm,g in named:
    got=None
    for m in range(2,7):
        if 2**m>len(g): break
        ok,p=descend(g,m)
        if ok: got=(m,p); break
    print("%-16s %s"%(nm, ("descends to 2^%d : %s"%(got[0],got[1])) if got else "*** NO DESCENT ***"))
    sys.stdout.flush()

print("\nrandom 3-connected cubic girth-5 graphs\n")
random.seed(5); bad=[]
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
        res=None
        for m in range(2,7):
            if 2**m>n: break
            ok,p=descend(adj,m)
            if ok: res=(m,p); break
        print("rand n=%-3d #%d  %s"%(n,got,("2^%d : %s"%(res[0],res[1])) if res else "*** NO DESCENT ***"))
        if res is None: bad.append((n,got))
        sys.stdout.flush()
print("\nfailures:",bad)
