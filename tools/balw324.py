"""balw324.py -- attack on (BAL-W).

(BAL-W) for SOME colouring with the three terminals monochromatic, W = B - N contains
        vertices of BOTH colours carrying fans to the three terminals inside W.
By 321.3 this implies (BAL), hence (N2).
Tested three ways: with the minimum-|N| monochromatic colouring; over a search of such
colourings; and with a note on whether W is even connected, since deleting N from a
2-connected graph can disconnect it and a component missing a terminal carries no fan.
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
def all_simple(a, s, t, banned):
    out=[]; p=[s]; on=set(banned)|{s}
    def rec(u):
        for w in a[u]:
            if w==t: out.append(list(p)+[t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out
def has_fan(a, v, terms):
    if v not in a or len(a[v])<3: return False
    for p0 in all_simple(a, v, terms[0], set(terms)-{terms[0]}):
        u0=set(p0[1:])
        for p1 in all_simple(a, v, terms[1], (u0-{v}) | {terms[2]}):
            u1=u0|set(p1[1:])
            if all_simple(a, v, terms[2], (u1-{v})): return True
    return False
def ncomp(a, vs):
    seen=set(); c=0
    for s in vs:
        if s in seen: continue
        c+=1; seen.add(s); q=deque([s])
        while q:
            u=q.popleft()
            for w in a.get(u,()):
                if w not in seen: seen.add(w); q.append(w)
    return c
def mono_colourings(a, terms, rng, k):
    out=[]; seen=set()
    for _ in range(k*10):
        col={v:rng.randint(0,1) for v in a}
        for x in terms: col[x]=0
        if rng.random()<0.7:
            imp=True
            while imp:
                imp=False
                for v in a:
                    if v in terms: continue
                    nc=sum(1 for w in a[v] if col[w]==col[v])
                    if nc>len(a[v])-nc: col[v]^=1; imp=True
        key=tuple(sorted((str(v),col[v]) for v in col))
        if key in seen: continue
        seen.add(key); out.append(dict(col))
        if len(out)>=k: break
    return out

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 103)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 26):
    a,terms=make_bridge(rng.choice([10,12]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    t=tau_free(E)
    if t<2: continue
    st['bridges tau>=2']+=1
    cols=mono_colourings(a,terms,rng,14)
    cols.sort(key=lambda c: sum(1 for e in E if c[tuple(e)[0]]==c[tuple(e)[1]]))
    okmin=None; okany=False; discon=0
    for i,col in enumerate(cols):
        N={e for e in E if col[tuple(e)[0]]==col[tuple(e)[1]]}
        if not N: continue
        aW=adj_of(E-N)
        if ncomp(aW, list(a))>1: discon+=1
        F={v for v in a if v not in terms and has_fan(aW,v,terms)}
        good=len({col[v] for v in F})==2
        if i==0: okmin=good
        if good: okany=True
    st['min-|N| colouring works' if okmin else 'min-|N| colouring FAILS']+=1
    st['(BAL-W) holds for SOME colouring' if okany else '(BAL-W) FAILS for all tried']+=1
    st['W disconnected for >=1 colouring']+=1 if discon else 0
for k in sorted(st): print('  %-40s %s'%(k,st[k]))
