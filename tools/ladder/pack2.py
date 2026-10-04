"""Push the packing measurement to larger n to test whether the supply of
pairwise-disjoint-span chords grows linearly (reaching 7 near n=40) or only
logarithmically.  Section 166.1 extrapolated from just the values 3 and 4."""
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

def build(n,pairs):
    adj=[[] for _ in range(n)]
    for i in range(n):
        adj[i].append((i+1)%n); adj[(i+1)%n].append(i)
    for a,b in pairs: adj[a].append(b); adj[b].append(a)
    return adj

def scan(adj, cap=400000):
    """iterate cycles; return best count of pairwise-disjoint-span chords"""
    n=len(adj); best=0; bestlen=0; cnt=0
    adjset=[set(a) for a in adj]
    for s in range(n):
        stack=[(s,[s],1<<s)]
        while stack:
            u,path,mask=stack.pop()
            for v in adj[u]:
                if v==s and len(path)>=3:
                    if path[1]<path[-1]:
                        cnt+=1
                        if cnt>cap: return best,bestlen,cnt
                        L=len(path); pos={w:i for i,w in enumerate(path)}
                        ce=set()
                        for i in range(L):
                            a,b=path[i],path[(i+1)%L]; ce.add((min(a,b),max(a,b)))
                        iv=[]
                        for w in path:
                            for z in adjset[w]:
                                if z in pos and w<z and (w,z) not in ce:
                                    p,q=pos[w],pos[z]
                                    if p>q: p,q=q,p
                                    inn=q-p; out=L-inn
                                    iv.append((p,q) if inn<=out else (q,p+L))
                        iv.sort(key=lambda t:t[1])
                        k=0; last=-1
                        for (l,r) in iv:
                            if l>=last: k+=1; last=r
                        if k>best: best=k; bestlen=L
                elif v>s and not (mask>>v)&1:
                    stack.append((v,path+[v],mask|(1<<v)))
    return best,bestlen,cnt

random.seed(7)
print("max pairwise-disjoint-span chords, maximised over all cycles")
print("n   : values (one per graph)          [166.1 predicted linear growth ~ n/5.5]\n")
for n in [24,26,28,30,32]:
    vals=[]; got=0; tries=0
    while got<3 and tries<60000:
        tries+=1
        perm=list(range(n)); random.shuffle(perm)
        pr=[(perm[2*i],perm[2*i+1]) for i in range(n//2)]
        if any((b-a)%n in (1,n-1) for a,b in pr): continue
        if len(set(tuple(sorted(p)) for p in pr))!=n//2: continue
        adj=build(n,pr)
        if any(len(x)!=3 for x in adj): continue
        if girth(adj)<5 or not three_conn(adj): continue
        got+=1
        b,bl,c=scan(adj)
        vals.append((b,bl,c))
    print("n=%-3d: %s      (best,cyclelen,#cycles)" % (n,vals))
    sys.stdout.flush()
