import itertools, random, sys
from collections import deque

def girth(adj):
    n=len(adj); best=10**9
    for s in range(n):
        dist={s:0}; par={s:None}; q=deque([s])
        while q:
            u=q.popleft()
            for w in adj[u]:
                if w not in dist:
                    dist[w]=dist[u]+1; par[w]=u; q.append(w)
                elif w!=par[u]:
                    best=min(best,dist[u]+dist[w]+1)
    return best

def named():
    G={}
    # Petersen
    adj=[[] for _ in range(10)]
    E=[(0,1),(1,2),(2,3),(3,4),(4,0),(5,7),(7,9),(9,6),(6,8),(8,5),(0,5),(1,6),(2,7),(3,8),(4,9)]
    G['Petersen']=(10,E)
    # Mobius-Kantor GP(8,3)
    E=[]
    for i in range(8):
        E.append((i,(i+1)%8)); E.append((8+i,8+(i+3)%8)); E.append((i,8+i))
    G['MobiusKantor']=(16,E)
    # Pappus
    E=[]
    pap=[[1,5,7],[2,10,14],[3,15,17],[4,8,12],[6,16,0],[9,13,11]]
    # use LCF for Pappus: [5,7,-7,7,-7,-5]^3 on 18
    lcf=[5,7,-7,7,-7,-5]*3; n=18; E=[(i,(i+1)%n) for i in range(n)]
    for i,s in enumerate(lcf):
        a,b=i,(i+s)%n
        if (min(a,b),max(a,b)) not in [(min(x,y),max(x,y)) for x,y in E]:
            E.append((a,b))
    G['Pappus']=(18,E)
    # Desargues LCF [5,-5,9,-9]^5 on 20
    lcf=[5,-5,9,-9]*5; n=20; E=[(i,(i+1)%n) for i in range(n)]
    for i,s in enumerate(lcf):
        a,b=i,(i+s)%n
        if (min(a,b),max(a,b)) not in [(min(x,y),max(x,y)) for x,y in E]:
            E.append((a,b))
    G['Desargues']=(20,E)
    # Dodecahedron LCF [10,7,4,-4,-7,10,-4,7,-7,4]^2 on 20
    lcf=[10,7,4,-4,-7,10,-4,7,-7,4]*2; n=20; E=[(i,(i+1)%n) for i in range(n)]
    for i,s in enumerate(lcf):
        a,b=i,(i+s)%n
        if (min(a,b),max(a,b)) not in [(min(x,y),max(x,y)) for x,y in E]:
            E.append((a,b))
    G['Dodecahedron']=(20,E)
    return G

def randcubic(n,rng,tries=4000):
    for _ in range(tries):
        stubs=[]
        for v in range(n): stubs += [v]*3
        rng.shuffle(stubs)
        E=set(); ok=True
        for i in range(0,len(stubs),2):
            a,b=stubs[i],stubs[i+1]
            if a==b or (min(a,b),max(a,b)) in E: ok=False; break
            E.add((min(a,b),max(a,b)))
        if not ok: continue
        adj=[[] for _ in range(n)]
        for a,b in E: adj[a].append(b); adj[b].append(a)
        if girth(adj)>=5: return sorted(E)
    return None

def build(n,E):
    adj=[[] for _ in range(n)]
    for a,b in E: adj[a].append(b); adj[b].append(a)
    return adj

def cycles(adj,n,maxlen):
    """all chordless cycles up to maxlen, as vertex lists"""
    out=[]
    for s in range(n):
        stack=[(s,[s],{s})]
        while stack:
            u,path,seen=stack.pop()
            for w in adj[u]:
                if w==s and len(path)>=5:
                    if min(path)==s and path[1]<path[-1]:
                        out.append(list(path))
                elif w not in seen and w>s and len(path)<maxlen:
                    stack.append((w,path+[w],seen|{w}))
    # keep chordless
    res=[]
    for cyc in out:
        L=len(cyc); pos={v:i for i,v in enumerate(cyc)}; ch=False
        for i,v in enumerate(cyc):
            for w in adj[v]:
                if w in pos:
                    d=abs(pos[v]-pos[w]); d=min(d,L-d)
                    if d>1: ch=True
        if not ch: res.append(cyc)
    return res

def allpathlens(adj,V,s,t,cap):
    """set of lengths of simple s-t paths inside vertex set V"""
    Vs=set(V); out=set()
    stack=[(s,{s},0)]
    while stack:
        u,seen,d=stack.pop()
        if d>cap: continue
        for w in adj[u]:
            if w not in Vs: continue
            if w==t: out.add(d+1)
            elif w not in seen:
                stack.append((w,seen|{w},d+1))
    return out

def main():
    rng=random.Random(20260805)
    graphs=[]
    for nm,(n,E) in named().items(): graphs.append((nm,n,E))
    for n in (14,16,18,20):
        for k in range(6):
            E=randcubic(n,rng)
            if E: graphs.append(("rand%d_%d"%(n,k),n,E))
    tot=0; cand=0; far=0
    for nm,n,E in graphs:
        adj=build(n,E)
        for C in cycles(adj,n,min(n,16)):
            L=len(C)
            if L<5: continue             # any chordless cycle: test the CORE statement
            Cs=set(C)
            rest=[v for v in range(n) if v not in Cs]
            if not rest: continue
            # components of G - C
            seen=set(); comps=[]
            for v in rest:
                if v in seen: continue
                q=deque([v]); comp=set([v]); seen.add(v)
                while q:
                    u=q.popleft()
                    for w in adj[u]:
                        if w in Cs or w in comp: continue
                        comp.add(w); seen.add(w); q.append(w)
                comps.append(comp)
            for comp in comps:
                terms=sorted(v for v in comp if any(w in Cs for w in adj[v]))
                if len(terms)<3: continue
                # distinct attachment neighbours (non-degenerate)
                att=[[w for w in adj[v] if w in Cs] for v in terms]
                if any(len(a)!=1 for a in att): continue
                if len({a[0] for a in att})!=len(terms): continue
                tot+=1
                cap=len(comp)
                lens={}
                ok=True
                mind=10**9
                for a,b in itertools.combinations(terms,2):
                    S=allpathlens(adj,comp,a,b,cap)
                    if not S: ok=False; break
                    mind=min(mind,min(S))
                    ev=sorted(x for x in S if x%2==0)
                    if len({x%4 for x in ev})>1: ok=False; break
                if not ok: continue
                cand+=1
                D=L//2-1
                tag="FAR" if mind>=D else "near"
                if mind>=D: far+=1
                print("CANDIDATE %-14s |C|=%d  s=%d  |B|=%d  mindist=%d  D=%d  %s"%(nm,L,len(terms),len(comp),mind,D,tag))
    print("\nbridges examined      :",tot)
    print("all-even-congruent    :",cand)
    print("AND terminals far apart:",far)

main()
