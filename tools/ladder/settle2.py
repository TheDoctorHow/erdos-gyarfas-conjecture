import itertools, random
from collections import deque

def girth(adj,n):
    best=10**9
    for s in range(n):
        dist={s:0}; par={s:None}; q=deque([s])
        while q:
            u=q.popleft()
            for w in adj[u]:
                if w not in dist:
                    dist[w]=dist[u]+1; par[w]=u; q.append(w)
                elif w!=par[u]: best=min(best,dist[u]+dist[w]+1)
    return best

def twoconn(adj,n):
    def conn(skip):
        st=[v for v in range(n) if v!=skip]
        if not st: return True
        seen={st[0]}; q=deque([st[0]])
        while q:
            u=q.popleft()
            for w in adj[u]:
                if w!=skip and w not in seen: seen.add(w); q.append(w)
        return len(seen)==n-(1 if skip is not None else 0)
    if not conn(None): return False
    return all(conn(v) for v in range(n))

def gen(n,s,rng,tries=300):
    """subcubic: s vertices of degree 2 (terminals 0..s-1), n-s of degree 3"""
    for _ in range(tries):
        stubs=[]
        for v in range(n): stubs += [v]*(2 if v<s else 3)
        if len(stubs)%2: return None
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
        if not twoconn(adj,n): continue
        return adj
    return None

def allpathlens(adj,s,t):
    out=set(); stack=[(s,{s},0)]
    while stack:
        u,seen,d=stack.pop()
        for w in adj[u]:
            if w==t: out.add(d+1)
            elif w not in seen: stack.append((w,seen|{w},d+1))
    return out

def main():
    rng=random.Random(7)
    tested=0; withint=0; hits=0; hits_int=0
    examples=[]
    for n in range(6,17):
        for s in range(3,min(n,9)+1):
            for trial in range(4000):
                adj=gen(n,s,rng,tries=1)
                if adj is None: continue
                tested+=1
                interior = n-s
                if interior>0: withint+=1
                ok=True
                for a,b in itertools.combinations(range(s),2):
                    S=allpathlens(adj,a,b)
                    if not S: ok=False; break
                    ev=[x for x in S if x%2==0]
                    if len({x%4 for x in ev})>1: ok=False; break
                if ok:
                    hits+=1
                    if interior>0:
                        hits_int+=1
                        if len(examples)<6:
                            examples.append((n,s,sorted((a,b) for a in range(n) for b in adj[a] if a<b)))
    print("subcubic girth>=5 2-connected bridges tested :",tested)
    print("  of which with >=1 INTERIOR vertex          :",withint)
    print("satisfying 'all even connecting paths = mod 4':",hits)
    print("  of which WITH an interior vertex           :",hits_int)
    for n,s,E in examples:
        print("   EXAMPLE n=%d s=%d  %s"%(n,s,E))

main()
