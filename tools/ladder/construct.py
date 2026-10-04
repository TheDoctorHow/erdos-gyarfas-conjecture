"""CONSTRUCT a counterexample to (H1).

Target: C = a 20-cycle (v_2(20)=2, not a power of two), t = 2.
Bridge B = {b, b'} joined by an edge, with FOUR attachments all at positions
= 0 mod 4, namely 0,4 -> b and 8,12 -> b'.  Then
    leaf ears  x0-b-x4      length 2, value 2 - 4 = -2,  v_2 = 1 < 2   shallow
    cross ears x-b-b'-y     length 3, value 3 - (multiple of 4),  odd, v_2 = 0 < 2
so B contributes nothing of v_2 >= 2.  Remaining 16 cycle vertices get 8 chords;
a chord of span s contributes 1 - s, which must have v_2 < 2, i.e. 1 - s odd or
= 2 mod 4, i.e. s even or s = 3 mod 4.   Simply: need v_2(1-s) < 2 for every chord,
i.e. s NOT = 1 mod 4.  Search chord matchings with girth >= 5 and 3-connectivity.
"""
import itertools, random, sys

def v2(d):
    if d==0: return 10**9
    k=0; d=abs(d)
    while d%2==0: d//=2; k+=1
    return k
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
def shifts_of(adj,cyc):
    VC=set(cyc); L=len(cyc); pos={w:i for i,w in enumerate(cyc)}
    ce=set()
    for i in range(L):
        a,b=cyc[i],cyc[(i+1)%L]; ce.add((min(a,b),max(a,b)))
    ears=[]
    for u in cyc:
        for v in adj[u]:
            if v in VC and u<v and (min(u,v),max(u,v)) not in ce: ears.append((u,v,1))
    for u in cyc:
        for w in adj[u]:
            if w in VC: continue
            st=[(w,[w],{w})]
            while st:
                x,path,seen=st.pop()
                for y in adj[x]:
                    if y in VC:
                        if y!=u and u<y: ears.append((u,y,len(path)+1))
                    elif y not in seen: st.append((y,path+[y],seen|{y}))
    out=[]
    for (u,v,ln) in ears:
        p,q=pos[u],pos[v]; a1=(q-p)%L; a2=L-a1
        out += [ln-a1, ln-a2]
    return out

L=20; b=20; bp=21; n=22
base=[(i,(i+1)%L) for i in range(L)]
base += [(0,b),(4,b),(b,bp),(8,bp),(12,bp)]
free=[p for p in range(L) if p not in (0,4,8,12)]   # 16 positions need chords
print("free positions:",free,len(free))

def ok_span(s):
    s=min(s,L-s)
    return 4<=s<=L-4                                  # girth>=5 needs span>=4

random.seed(2026)
best=None
for trial in range(300000):
    fr=free[:]; random.shuffle(fr)
    pairs=[(fr[2*i],fr[2*i+1]) for i in range(len(fr)//2)]
    if any(not ok_span(abs(x-y)) for x,y in pairs): continue
    E=base+[tuple(sorted(p)) for p in pairs]
    if len(set(E))!=len(E): continue
    adj=[[] for _ in range(n)]
    for x,y in E: adj[x].append(y); adj[y].append(x)
    if any(len(a)!=3 for a in adj): continue
    if girth(adj)<5: continue
    if not three_conn(adj): continue
    cyc=list(range(L))
    ds=shifts_of(adj,cyc)
    t=v2(L)
    if max(v2(d) for d in ds)<t:
        best=(pairs,adj,ds); break
if best:
    pairs,adj,ds=best
    print("\n*** COUNTEREXAMPLE TO (H1) FOUND ***")
    print("n=%d, C = 0..19 (|C|=20, t=2), bridge {20,21} joined by an edge"%n)
    print("attachments: 0,4 -> 20 ; 8,12 -> 21")
    print("chords:",sorted(pairs))
    print("all shifts:",sorted(set(ds)))
    print("max v_2 among shifts =",max(v2(d) for d in ds),"  (t = 2) -> SHALLOW")
    print("girth =",girth(adj),"  3-connected =",three_conn(adj))
else:
    print("\nno shallow completion found in 300000 trials")
