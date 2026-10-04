"""exch326.py -- the extremal principle for (CONN+): maximise CONNECTIVITY, not the cut.

At a maximum cut every flip strictly decreases the cut, so no flip that repairs a 2-cut can
preserve maximality: the maximum cut is the wrong extremal choice (as 322.2 already hinted).
Instead do local search on the objective
     f(c) = number of branch vertices of W carrying a fan to the three terminals,
which by Menger is the number with 3 internally disjoint paths to the apex T*.
Question: does local search on f always reach a colouring witnessing (BAL-W)?
Also recorded: the structural fact that any 2-cut of W+ must be crossed by an N-edge.
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

def vertex_conn_ge3(a, s, t):
    """>=3 internally disjoint s-t paths, via unit vertex capacities and BFS augmenting"""
    cap=collections.defaultdict(int); adjf=collections.defaultdict(set)
    def add(u,v,c):
        adjf[u].add(v); adjf[v].add(u); cap[(u,v)]+=c
    for v in a:
        c = 10**6 if v in (s,t) else 1
        add((v,'i'),(v,'o'),c)
    for u in a:
        for w in a[u]:
            add((u,'o'),(w,'i'),10**6)
    src=(s,'o'); snk=(t,'i'); flow=0
    while flow<3:
        par={src:None}; q=deque([src]); found=False
        while q and not found:
            u=q.popleft()
            for w in adjf[u]:
                if w in par or cap[(u,w)]<=0: continue
                par[w]=u; 
                if w==snk: found=True; break
                q.append(w)
        if not found: return False
        z=snk
        while par[z] is not None:
            cap[(par[z],z)]-=1; cap[(z,par[z])]+=1; z=par[z]
        flow+=1
    return True

def fancount(a, terms, E, N, col):
    Wp=set(E-N)
    for x in terms: Wp.add(frozenset(('T*',x)))
    aW=adj_of(Wp)
    F=[v for v in a if v not in terms and v in aW and len(aW[v])>=3
       and vertex_conn_ge3(aW, v, 'T*')]
    return F

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 151)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 24):
    a,terms=make_bridge(rng.choice([10,12]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    if tau_free(E)<2: continue
    st['bridges tau>=2']+=1
    col={v:rng.randint(0,1) for v in a}
    for x in terms: col[x]=0
    def NN(c): return {e for e in E if c[tuple(e)[0]]==c[tuple(e)[1]]}
    best=len(fancount(a,terms,E,NN(col),col))
    imp=True; steps=0
    while imp and steps<60:
        imp=False; steps+=1
        for v in a:
            if v in terms: continue
            col[v]^=1
            val=len(fancount(a,terms,E,NN(col),col))
            if val>best: best=val; imp=True
            else: col[v]^=1
    F=fancount(a,terms,E,NN(col),col)
    both=len({col[v] for v in F})==2
    st['local search on connectivity reaches (BAL-W)' if both
       else '*** local search STUCK ***']+=1
    st['  fan-carrying vertices at optimum: %d'%len(F)]+=1
for k in sorted(st): print('  %-46s %s'%(k,st[k]))
