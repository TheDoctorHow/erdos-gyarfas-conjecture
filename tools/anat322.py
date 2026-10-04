"""anat322.py -- anatomy of the misses, with EXHAUSTIVE fan enumeration.

The misses of 318/319 were found with a randomized 160-trial fan search.  317 already
caught one confounded measurement; before drawing structural conclusions, enumerate ALL
fans exhaustively (n is only 13-17) and see whether the misses survive.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/yparity315.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def tau_free(E, cap=4):
    El=sorted(E,key=lambda s:sorted(map(str,s)))
    for k in range(cap+1):
        for S in itertools.combinations(El,k):
            if bipn(adj_of(set(El)-set(S))) is not None: return k
    return cap+1

def all_simple(a, s, t, banned):
    out=[]; p=[s]; on=set(banned)|{s}
    def rec(u):
        for w in a[u]:
            if w==t: out.append(list(p)+[t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out

def all_fans_sigma(a, v, terms):
    """exhaustive: set of leg-parity triples over ALL fans at v"""
    sig=set()
    P0=all_simple(a, v, terms[0], set(terms)-{terms[0]})
    for p0 in P0:
        u0=set(p0[1:])
        P1=all_simple(a, v, terms[1], (u0-{v}) | ({terms[2]}))
        for p1 in P1:
            u1=u0|set(p1[1:])
            P2=all_simple(a, v, terms[2], (u1-{v}))
            for p2 in P2:
                sig.add((( len(p0)-1)%2, (len(p1)-1)%2, (len(p2)-1)%2))
        if len(sig)==8: break
    return sig

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 77)
st=collections.Counter(); dumped=0
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 40):
    a,terms=make_bridge(rng.choice([10,12]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    t=tau_free(E)
    if t<2: continue
    best=None
    for trial in range(250):
        col={v:rng.randint(0,1) for v in a}
        for x in terms: col[x]=0
        imp=True
        while imp:
            imp=False
            for v in a:
                if v in terms: continue
                nc=sum(1 for w in a[v] if col[w]==col[v])
                if nc>len(a[v])-nc: col[v]^=1; imp=True
        N={e for e in E if col[tuple(e)[0]]==col[tuple(e)[1]]}
        if best is None or len(N)<len(best[1]): best=(dict(col),set(N))
    col,N=best; t3=len(N); d=t3-t
    # EXHAUSTIVE: sigma over all balanced fans
    SIG=set(); bal=[]
    for v in a:
        if v in terms: continue
        S=all_fans_sigma(a,v,terms)
        for tr in S:
            if len(set(tr))==1:
                SIG.add(tr[0]); bal.append((v,tr[0]))
    fires=(len(SIG)==2)
    st['bridges tau>=2']+=1
    st['delta=%d : %s'%(d,'FIRES' if fires else 'MISS ')]+=1
    if not fires and dumped<3:
        dumped+=1
        print('  --- MISS: n=%d tau=%d tau_3=%d delta=%d  sigma=%s'%(len(a),t,t3,d,sorted(SIG)))
        print('      balanced fans: %d, at vertices %s'
              %(len(bal), sorted({str(v) for v,_ in bal})))
        print('      colour split of balanced-fan centres: %s'
              %collections.Counter(col[v] for v,_ in bal))
        nofan=[v for v in a if v not in terms and not all_fans_sigma(a,v,terms)]
        print('      vertices carrying NO fan at all: %d of %d'%(len(nofan),len(a)-3))
for k in sorted(st): print('  %-26s %s'%(k,st[k]))
