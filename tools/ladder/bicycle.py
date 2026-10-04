import itertools, random
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
    return best

def twoconn(adj,n):
    def conn(skip):
        st=[v for v in range(n) if v!=skip]
        seen={st[0]}; q=deque([st[0]])
        while q:
            u=q.popleft()
            for w in adj[u]:
                if w!=skip and w not in seen: seen.add(w); q.append(w)
        return len(seen)==len(st)
    return all(conn(v) for v in range(n)) and conn(None)

def gen(n,s,rng):
    stubs=[]
    for v in range(n): stubs += [v]*(2 if v<s else 3)
    if len(stubs)%2: return None
    rng.shuffle(stubs)
    E=set()
    for i in range(0,len(stubs),2):
        a,b=stubs[i],stubs[i+1]
        if a==b or (min(a,b),max(a,b)) in E: return None
        E.add((min(a,b),max(a,b)))
    adj=[[] for _ in range(n)]
    for a,b in E: adj[a].append(b); adj[b].append(a)
    if girth(adj,n)<5: return None
    if not twoconn(adj,n): return None
    return sorted(E),adj

def rowreduce(rows,ncols):
    rows=[r for r in rows if r]; basis=[]
    for r in rows:
        for b in basis:
            low=b & -b
            if r & low: r ^= b
        if r: basis.append(r); basis.sort(key=lambda x:-(x & -x))
    # simple gaussian
    basis=[]
    for r in rows:
        cur=r
        for b in basis:
            cur=min(cur,cur^b) if False else cur
        basis.append(cur)
    return None

def gauss(vectors):
    """return a basis (list of ints) of the span, over F2"""
    piv={}
    for v in vectors:
        cur=v
        while cur:
            h=cur.bit_length()-1
            if h in piv: cur^=piv[h]
            else: piv[h]=cur; break
    return list(piv.values())

def spaces(E,n):
    m=len(E)
    # cut space: rows = for each vertex, the set of incident edges
    rows=[]
    for v in range(n):
        r=0
        for i,(a,b) in enumerate(E):
            if a==v or b==v: r|=1<<i
        rows.append(r)
    cut=gauss(rows)
    # cycle space = orthogonal complement of cut space (over F2, wrt standard form)
    # build by solving: x with incidence*x = 0
    # use gaussian elimination on the incidence matrix (rows=vertices, cols=edges)
    piv={}; mat=[]
    for r in rows:
        cur=r
        for h in sorted(piv,reverse=True):
            if cur>>h & 1: cur^=piv[h]
        if cur: piv[cur.bit_length()-1]=cur
    pivcols=set(piv.keys())
    free=[i for i in range(m) if i not in pivcols]
    cyc=[]
    for f in free:
        x=1<<f
        # back-substitute
        for h in sorted(piv):
            row=piv[h]
            if bin(row & x).count('1')%2: x^=1<<h
        cyc.append(x)
    return cut,cyc

def dim_inter(A,B):
    """dim of intersection of spans of A and B over F2"""
    a=gauss(A); b=gauss(B)
    ab=gauss(a+b)
    return len(a)+len(b)-len(ab)

def main():
    rng=random.Random(11)
    worst={}
    viol=0; tested=0
    for n in range(6,17):
        for s in range(3,min(n,9)+1):
            for _ in range(3000):
                g=gen(n,s,rng)
                if not g: continue
                E,adj=g
                interior=n-s
                if interior<1: continue
                tested+=1
                mu=len(E)-n+1
                cut,cyc=spaces(E,n)
                d=dim_inter(cut,cyc)
                key=(interior,)
                if key not in worst or d-(mu-2)>worst[key][0]:
                    worst[key]=(d-(mu-2),n,s,mu,d)
                if d>=mu-2:
                    viol+=1
                    if viol<=5:
                        print("  d >= mu-2 !  n=%d s=%d interior=%d mu=%d d=%d"%(n,s,interior,mu,d))
    print("bridges with >=1 interior vertex tested:",tested)
    print("with bicycle dimension d >= mu-2      :",viol)
    print()
    for k in sorted(worst):
        gap,n,s,mu,d=worst[k]
        print("  interior=%-3d  best case: n=%d s=%d mu=%d d=%d   (d-(mu-2)=%d)"%(k[0],n,s,mu,d,gap))

main()
