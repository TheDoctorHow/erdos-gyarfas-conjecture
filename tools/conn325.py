"""conn325.py -- (CONN) restated with an apex, and what a maximum cut gives.

A fan at v to b_1,b_2,b_3 inside W is exactly three internally disjoint v-t* paths in
     W+ = W together with a new apex t* joined to b_1, b_2, b_3,
so by MENGER it exists iff no 2-cut of W+ separates v from t*.  Suppressing the degree-2
vertices of W+ and asking for 3-connectivity is therefore the correct form of (CONN); the
version in 322.4 kept the terminals unsuppressed and so could never hold.
Measured: does some monochromatic colouring make the suppressed W+ 3-connected?
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split(chr(34)*3,2)[2])

def bipn(a):
    col={}
    for s in a:
        if s in col: continue
        col[s]=0;q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in col: col[w]=1-col[u];q.append(w)
                elif col[w]==col[u]: return None
    return col
def tau_free(E, cap=4):
    El=sorted(E,key=lambda s:sorted(map(str,s)))
    for k in range(cap+1):
        for S in itertools.combinations(El,k):
            if bipn(adj_of(set(El)-set(S))) is not None: return k
    return cap+1
def suppress(a):
    a={v:set(ws) for v,ws in a.items()}
    ch=True
    while ch:
        ch=False
        for v in list(a):
            if len(a.get(v,()))==2:
                x,y=tuple(a[v])
                if x==y or y in a[x]: continue
                a[x].discard(v); a[y].discard(v); a[x].add(y); a[y].add(x); del a[v]; ch=True
    return a
def kconn(a,k):
    vs=list(a)
    if len(vs)<=k: return False
    for S in itertools.combinations(vs,k-1):
        rem=[v for v in vs if v not in S]
        seen={rem[0]}; q=deque([rem[0]])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w in S or w in seen: continue
                seen.add(w); q.append(w)
        if len(seen)!=len(rem): return False
    return True

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 131)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 26):
    a,terms=make_bridge(rng.choice([10,12]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    if tau_free(E)<2: continue
    st['bridges tau>=2']+=1
    ok=False; maxcut_ok=False; mincross=[]
    for trial in range(60):
        col={v:rng.randint(0,1) for v in a}
        for x in terms: col[x]=0
        imp=True
        while imp:                      # local max cut, non-terminals only
            imp=False
            for v in a:
                if v in terms: continue
                nc=sum(1 for w in a[v] if col[w]==col[v])
                if nc>len(a[v])-nc: col[v]^=1; imp=True
        N={e for e in E if col[tuple(e)[0]]==col[tuple(e)[1]]}
        if not N: continue
        Wp=set(E-N)
        for x in terms: Wp.add(frozenset(('T*',x)))
        aW=adj_of(Wp)
        good=kconn(suppress(aW),3)
        if trial==0: maxcut_ok=good
        mincross.append(min(sum(1 for w in a[v] if col[w]!=col[v])
                            for v in a if v not in terms))
        if good: ok=True; break
    st['(CONN+) holds for some colouring' if ok else '(CONN+) FAILS for all tried']+=1
    st['a locally-maximum cut suffices' if maxcut_ok else 'local max cut not enough']+=1
    if mincross:
        st['min crossing-degree at non-terminals = %d'%min(mincross)]+=1
for k in sorted(st): print('  %-42s %s'%(k,st[k]))
