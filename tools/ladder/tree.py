import itertools
from collections import deque

def prufer_to_tree(seq,n):
    deg=[1]*n
    for x in seq: deg[x]+=1
    E=[]; s=list(seq)
    import heapq
    leaves=[i for i in range(n) if deg[i]==1]; heapq.heapify(leaves)
    ptr=0
    for x in s:
        leaf=heapq.heappop(leaves)
        E.append((leaf,x)); deg[x]-=1
        if deg[x]==1: heapq.heappush(leaves,x)
    u=heapq.heappop(leaves); v=heapq.heappop(leaves); E.append((u,v))
    return E

def steiner_trees(k):
    """all trees on k leaves (0..k-1, degree 1) + (k-2) internal (degree 3)"""
    n=2*k-2; internals=list(range(k,n))
    out=[]
    for seq in set(itertools.permutations(sum([[i,i] for i in internals],[]))):
        E=prufer_to_tree(list(seq),n)
        deg=[0]*n
        for a,b in E: deg[a]+=1; deg[b]+=1
        if all(deg[i]==1 for i in range(k)) and all(deg[i]==3 for i in internals):
            out.append(E)
    # dedup
    seen=set(); res=[]
    for E in out:
        key=tuple(sorted((min(a,b),max(a,b)) for a,b in E))
        if key not in seen: seen.add(key); res.append(E)
    return res

def pathedges(E,n,a,b):
    adj={i:[] for i in range(n)}
    for idx,(x,y) in enumerate(E): adj[x].append((y,idx)); adj[y].append((x,idx))
    par={a:(None,None)}; q=deque([a])
    while q:
        u=q.popleft()
        for w,idx in adj[u]:
            if w not in par: par[w]=(u,idx); q.append(w)
    out=[]; cur=b
    while cur!=a:
        p,idx=par[cur]; out.append(idx); cur=p
    return out

def test(k,positions):
    n=2*k-2
    trees=steiner_trees(k)
    ne=n-1
    sat=[]
    for E in trees:
        pe={}
        for a,b in itertools.combinations(range(k),2):
            pe[(a,b)]=pathedges(E,n,a,b)
        for assign in itertools.product(range(4),repeat=ne):
            ok=True
            for (a,b),eidx in pe.items():
                d=sum(assign[i] for i in eidx)%4
                if d==(positions[b]-positions[a])%4: ok=False; break
            if ok:
                sat.append((E,assign)); break
    return len(trees),sat

for k in (4,5,6):
    pos=list(range(k))
    nt,sat=test(k,pos)
    print("k=%d marked points, consecutive positions %s : %d Steiner topologies, %d satisfiable"
          %(k,pos,nt,len(sat)))
    if sat and k<=5:
        E,a=sat[0]
        print("    e.g. tree",sorted(E)," residues",a)

# also: non-consecutive position patterns for k=5
print()
best=None
for pos in itertools.combinations(range(0,9),5):
    nt,sat=test(5,list(pos))
    if len(sat)==0:
        print("UNSATISFIABLE for positions",pos)
        best=pos; break
if best is None:
    print("every 5-position pattern in range 0..8 is satisfiable for some tree")
