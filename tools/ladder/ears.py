"""EAR-SUBSET TEST.

LEMMA (elementary, proved). Let C be a cycle and P_1..P_t ears (paths with both
ends on C, interiors disjoint from C and from each other) whose spanned arcs
A_1..A_t are PAIRWISE DISJOINT sub-arcs of C.  Put delta_j = |P_j| - |A_j|.
Then for every S subset of [t], G has a cycle of length |C| + sum_{j in S} delta_j.
(Replace A_j by P_j for j in S; disjointness of the arcs keeps it a cycle.)

So G has a cycle of length 0 mod 8 as soon as some subset sum of the delta_j
equals -|C| mod 8.

delta_j is ODD  <=>  P_j together with A_j is an ODD cycle.

TEST: for Hamiltonian cubic graphs G = C_n + perfect matching, the ears are the
chords (|P_j| = 1), so delta_j = 1 - |A_j|.  We take a maximum set of chords with
pairwise disjoint spans (interval scheduling on the shorter arc) and ask whether
the subset sums of their deltas hit -n mod 8, and with how few chords.
"""
import random, itertools

def build(n, pairs):
    adj=[[] for _ in range(n)]
    for i in range(n):
        adj[i].append((i+1)%n); adj[(i+1)%n].append(i)
    for a,b in pairs:
        adj[a].append(b); adj[b].append(a)
    return adj

def girth(adj):
    n=len(adj); best=10**9
    for s in range(n):
        dist=[-1]*n; par=[-1]*n; dist[s]=0; q=[s]
        for u in q:
            for v in adj[u]:
                if dist[v]<0: dist[v]=dist[u]+1; par[v]=u; q.append(v)
                elif v!=par[u]: best=min(best,dist[u]+dist[v]+1)
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
    n=len(adj)
    return all(conn_after(adj,set(p)) for p in itertools.combinations(range(n),2))

def analyse(n, pairs, name):
    # chords as intervals on the Hamiltonian cycle 0..n-1
    ivals=[]
    for a,b in pairs:
        x,y=min(a,b),max(a,b)
        inner=y-x            # arc going x -> y forward
        outer=n-inner
        if inner<=outer: ivals.append((x,y,inner))
        else:            ivals.append((y,x+n,outer))   # wrap
    # interval scheduling: maximum set of pairwise disjoint arcs (by right endpoint)
    ivals.sort(key=lambda t:t[1])
    chosen=[]; last=-1
    for (l,r,ln) in ivals:
        if l>=last:
            chosen.append((l,r,ln)); last=r
    deltas=[(1-ln)%8 for (l,r,ln) in chosen]
    target=(-n)%8
    # minimum subset size hitting target
    best=None
    reach={0:0}
    for d in deltas:
        new=dict(reach)
        for s,c in reach.items():
            t=(s+d)%8
            if t not in new or new[t]>c+1: new[t]=c+1
        reach=new
    hit = target in reach
    best = reach.get(target)
    return dict(name=name, n=n, girth=girth(build(n,pairs)), nchords=len(pairs),
                disjoint=len(chosen), target=target, deltas=sorted(set(deltas)),
                hit=hit, minsubset=best, hamlen_mod8=n%8)

print("EAR-SUBSET TEST on 3-connected cubic girth->=5 Hamiltonian graphs")
print("target = (-n) mod 8 ; 'hit' = some subset of disjoint-chord deltas reaches it\n")
random.seed(20260805)
fails=0; tot=0; maxneed=0
for n in range(12,41,2):
    got=0; tries=0
    while got<6 and tries<40000:
        tries+=1
        perm=list(range(n)); random.shuffle(perm)
        pr=[(perm[2*i],perm[2*i+1]) for i in range(n//2)]
        if any((b-a)%n in (1,n-1) for a,b in pr): continue
        adj=build(n,pr)
        if any(len(x)!=3 for x in adj): continue
        if len(set(tuple(sorted(p)) for p in pr))!=n//2: continue
        if girth(adj)<5: continue
        if not three_conn(adj): continue
        got+=1; tot+=1
        r=analyse(n,pr,"rand n=%d #%d"%(n,got))
        if not r["hit"]: fails+=1
        else: maxneed=max(maxneed,r["minsubset"])
        if got<=2 or not r["hit"]:
            print("%-14s n=%-3d g=%d chords=%-3d disjoint-span chords=%-3d "
                  "target=%d deltas=%s  HIT=%-5s minsubset=%s"
                  % (r["name"],r["n"],r["girth"],r["nchords"],r["disjoint"],
                     r["target"],r["deltas"],r["hit"],r["minsubset"]))
print()
print("graphs tested: %d   FAILURES (no subset reaches target): %d" % (tot,fails))
print("largest minimum subset size needed:", maxneed)
