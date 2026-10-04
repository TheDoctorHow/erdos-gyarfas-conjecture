import random, itertools, sys
from collections import deque

def girth(adj,n):
    best=10**9
    for s in range(n):
        dist={s:0}; par={s:None}; q=deque([s])
        while q:
            u=q.popleft()
            for w in adj[u]:
                if w not in dist: dist[w]=dist[u]+1; par[w]=u; q.append(w)
                elif w!=par[u]: best=min(best,dist[u]+dist[w]+1)
        if best<=4: return best
    return best

def randcubic(n,rng,tries=800):
    for _ in range(tries):
        stubs=[]
        for v in range(n): stubs+=[v]*3
        rng.shuffle(stubs)
        E=set(); ok=True
        for i in range(0,len(stubs),2):
            a,b=stubs[i],stubs[i+1]
            if a==b or (min(a,b),max(a,b)) in E: ok=False; break
            E.add((min(a,b),max(a,b)))
        if not ok: continue
        adj=[[] for _ in range(n)]
        for a,b in E: adj[a].append(b); adj[b].append(a)
        if girth(adj,n)<5: continue
        return sorted(E),adj
    return None

def dist(adj,n,s):
    d=[-1]*n; d[s]=0; q=deque([s])
    while q:
        u=q.popleft()
        for w in adj[u]:
            if d[w]<0: d[w]=d[u]+1; q.append(w)
    return d

def cycles_upto(adj,n,L):
    out=set()
    for s in range(n):
        stack=[(s,(s,),{s})]
        while stack:
            u,path,seen=stack.pop()
            for w in adj[u]:
                if w==s and len(path)>=4:
                    out.add(frozenset(path))
                elif w not in seen and w>s and len(path)<L:
                    stack.append((w,path+(w,),seen|{w}))
    return out

def cycedges(adj,C):
    # edges of the induced cycle on vertex set C (2-regular in cubic graph)
    E=[]
    for u in C:
        for w in adj[u]:
            if w in C and u<w: E.append((u,w))
    return E

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 4)
    tot=0; withnu3=0; less=[]
    for nH in (18,20,22,24,26):
        for trial in range(120):
            g=randcubic(nH,rng)
            if not g: continue
            E,adj=g
            D=[dist(adj,nH,v) for v in range(nH)]
            def edist(e,f):
                return min(D[a][b] for a in e for b in f)
            # matchings of size 3 with pairwise edge-distance >= 3
            cand=[]
            for trip in itertools.combinations(E,3):
                if all(edist(trip[i],trip[j])>=3 for i in range(3) for j in range(i+1,3)):
                    cand.append(trip)
            if not cand: continue
            rng.shuffle(cand)
            for M in cand[:6]:
                Ms=set(M)
                tot+=1
                # negative cycles: |Z| + |Z cap M| odd
                neg=[]
                for C in cycles_upto(adj,nH,14):
                    ce=cycedges(adj,C)
                    if len(ce)!=len(C): continue
                    k=sum(1 for e in ce if e in Ms)
                    if (len(C)+k)%2==1: neg.append(frozenset(C))
                # max packing of disjoint negative cycles (cap at 3)
                best=0
                def rec(i,used,cnt):
                    nonlocal best
                    if cnt>best: best=cnt
                    if best>=3 or i>=len(neg): return
                    for j in range(i,len(neg)):
                        if not (neg[j]&used): rec(j+1,used|neg[j],cnt+1)
                rec(0,frozenset(),0)
                if best>=3: withnu3+=1
                else: less.append((nH,best,len(neg)))
    print("valid surviving configurations tested :",tot)
    print("  with 3 disjoint negative cycles     :",withnu3)
    print("  with fewer than 3                   :",len(less))
    from collections import Counter
    if less:
        print("  breakdown (nH, nu_odd, #neg cycles found):")
        for k,v in sorted(Counter(less).items())[:20]: print("    ",k,"x",v)

main()
