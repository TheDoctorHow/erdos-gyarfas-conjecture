"""Direct test of the THREE-RESIDUE EAR STATEMENT (167.3).

For every cycle C of a cubic graph, list its chords; for a chord spanning arc A
put delta = 1 - |A| (mod 8) and v = v_2(delta) (v=3 means delta = 0 mod 8).
Ask: is there a cycle C carrying chords with pairwise DISJOINT arcs whose
valuations include 0,1,2  (general case)  or 1,2 (bipartite case)?
Arcs are allowed to be either of the two arcs of the chord.
"""
import random, itertools, sys

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

def bipartite(adj):
    n=len(adj); col=[-1]*n
    for s in range(n):
        if col[s]>=0: continue
        col[s]=0; q=[s]
        for u in q:
            for v in adj[u]:
                if col[v]<0: col[v]=1-col[u]; q.append(v)
                elif col[v]==col[u]: return False
    return True

def build(n,pairs):
    adj=[[] for _ in range(n)]
    for i in range(n):
        adj[i].append((i+1)%n); adj[(i+1)%n].append(i)
    for a,b in pairs: adj[a].append(b); adj[b].append(a)
    return adj

def v2(d):
    d%=8
    if d==0: return 3
    k=0
    while d%2==0: d//=2; k+=1
    return k

def disjoint(a,b):
    # arcs given as (start,end) half-open in [0,2L), disjoint on the circle of length L
    (s1,e1),(s2,e2)=a,b
    return e1<=s2 or e2<=s1

def test_graph(adj, need):
    """need = set of valuations required, e.g. {0,1,2} or {1,2}"""
    n=len(adj); adjset=[set(x) for x in adj]
    for s in range(n):
        stack=[(s,[s],1<<s)]
        while stack:
            u,path,mask=stack.pop()
            for v in adj[u]:
                if v==s and len(path)>=3:
                    if path[1]<path[-1]:
                        L=len(path); pos={w:i for i,w in enumerate(path)}
                        ce=set()
                        for i in range(L):
                            a,b=path[i],path[(i+1)%L]; ce.add((min(a,b),max(a,b)))
                        opts=[]   # (arc as interval, valuation)
                        for w in path:
                            for z in adjset[w]:
                                if z in pos and w<z and (w,z) not in ce:
                                    p,q=pos[w],pos[z]
                                    if p>q: p,q=q,p
                                    for (st,en,ln) in ((p,q,q-p),(q,p+L,L-(q-p))):
                                        opts.append(((st,en), v2(1-ln)))
                        # brute force over small subsets
                        k=len(need)
                        for combo in itertools.combinations(range(len(opts)),k):
                            vals={opts[i][1] for i in combo}
                            if vals!=need: continue
                            ok=True
                            for i,j in itertools.combinations(combo,2):
                                if not disjoint(opts[i][0],opts[j][0]): ok=False;break
                            if ok: return True, L
                elif v>s and not (mask>>v)&1:
                    stack.append((v,path+[v],mask|(1<<v)))
    return False, None

random.seed(11)
print("TEST of 167.3 -- three-residue (general) / two-residue (bipartite) ear statement\n")
okA=0; totA=0; failA=[]
for n in range(12,25,2):
    got=0; tries=0
    while got<4 and tries<60000:
        tries+=1
        perm=list(range(n)); random.shuffle(perm)
        pr=[(perm[2*i],perm[2*i+1]) for i in range(n//2)]
        if any((b-a)%n in (1,n-1) for a,b in pr): continue
        if len(set(tuple(sorted(p)) for p in pr))!=n//2: continue
        adj=build(n,pr)
        if any(len(x)!=3 for x in adj): continue
        if girth(adj)<5 or not three_conn(adj): continue
        got+=1; totA+=1
        bip=bipartite(adj)
        need={1,2} if bip else {0,1,2}
        ok,L=test_graph(adj,need)
        if ok: okA+=1
        else: failA.append((n,bip))
        print("n=%-3d bip=%-5s need=%s -> %s%s" % (n,bip,sorted(need),ok,
              ("  (witness cycle length %d)"%L) if ok else ""))
        sys.stdout.flush()
print()
print("graphs: %d   satisfying the statement: %d   failures: %s" % (totA,okA,failA))
