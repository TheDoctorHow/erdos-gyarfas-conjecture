"""Measure the packing constant: over ALL cycles C of a cubic girth->=5 graph,
the maximum number of chords of C with PAIRWISE DISJOINT spans.
Section 165.3 needs about 8.  Also record whether the ear-subset criterion
(some subset of deltas = -|C| mod 8) is met for the BEST cycle."""
import random, itertools

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

def all_cycles(adj):
    n=len(adj); out=[]
    for s in range(n):
        stack=[(s,[s],1<<s)]
        while stack:
            u,path,mask=stack.pop()
            for v in adj[u]:
                if v==s and len(path)>=3:
                    if path[1]<path[-1]: out.append(path[:])
                elif v>s and not (mask>>v)&1:
                    stack.append((v,path+[v],mask|(1<<v)))
    return out

def max_disjoint_chords(adj, cyc):
    L=len(cyc); pos={v:i for i,v in enumerate(cyc)}
    onc=set(cyc); cyc_edges=set()
    for i in range(L):
        a,b=cyc[i],cyc[(i+1)%L]; cyc_edges.add((min(a,b),max(a,b)))
    ivals=[]
    for i,u in enumerate(cyc):
        for v in adj[u]:
            if v in onc and (min(u,v),max(u,v)) not in cyc_edges and u<v:
                p,q=pos[u],pos[v]
                if p>q: p,q=q,p
                inner=q-p; outer=L-inner
                if inner<=outer: ivals.append((p,q,inner))
                else: ivals.append((q,p+L,outer))
    ivals.sort(key=lambda t:t[1])
    chosen=[]; last=-1
    for (l,r,ln) in ivals:
        if l>=last: chosen.append((l,r,ln)); last=r
    return len(chosen), [(1-ln)%8 for (l,r,ln) in chosen]

def report(name, adj):
    n=len(adj); best=0; bestinfo=None; anyhit=False
    for cyc in all_cycles(adj):
        k,deltas=max_disjoint_chords(adj,cyc)
        if k>best:
            best=k; bestinfo=(len(cyc),deltas)
        # criterion check
        reach={0}
        for d in deltas:
            reach |= {(s+d)%8 for s in reach}
        if (-len(cyc))%8 in reach: anyhit=True
    print("%-14s n=%-3d  max disjoint-span chords over all cycles = %-2d  "
          "(best cycle len %s, deltas %s)   some cycle meets criterion: %s"
          % (name,n,best,bestinfo[0] if bestinfo else "-",
             sorted(set(bestinfo[1])) if bestinfo else "-", anyhit))
    return best, anyhit

def build(n,pairs):
    adj=[[] for _ in range(n)]
    for i in range(n):
        adj[i].append((i+1)%n); adj[(i+1)%n].append(i)
    for a,b in pairs: adj[a].append(b); adj[b].append(a)
    return adj

random.seed(20260805)
print("PACKING MEASUREMENT (section 165.3 needs ~8 disjoint-span ears)\n")
allbest=[]; allhit=0; tot=0
for n in range(12,23,2):
    got=0; tries=0
    while got<4 and tries<40000:
        tries+=1
        perm=list(range(n)); random.shuffle(perm)
        pr=[(perm[2*i],perm[2*i+1]) for i in range(n//2)]
        if any((b-a)%n in (1,n-1) for a,b in pr): continue
        if len(set(tuple(sorted(p)) for p in pr))!=n//2: continue
        adj=build(n,pr)
        if any(len(x)!=3 for x in adj): continue
        if girth(adj)<5 or not three_conn(adj): continue
        got+=1; tot+=1
        b,h=report("rand n=%d #%d"%(n,got),adj)
        allbest.append((n,b))
        if h: allhit+=1
print()
print("graphs:",tot,"  cycles meeting the ear-subset criterion in at least one cycle:",allhit)
print("max disjoint-span chords by n:", allbest)
