"""even328.py -- what characterises a vertex carrying an EVEN balanced fan?

Target (325.4): some vertex of B carries a fan whose three legs are all even.
Rather than guess a characterisation, test several at once against exhaustive enumeration:
  H1  d(v,b_m) even for all three m
  H2  both parities of v-b_m path are achievable for all three m  (parity-flexible)
  H3  v lies on some odd cycle
and record how the sigma=0 set relates to each.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split(chr(34)*3,2)[2])

def all_simple(a, s, t, banned):
    out=[]; p=[s]; on=set(banned)|{s}
    def rec(u):
        for w in a[u]:
            if w==t: out.append(list(p)+[t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out
def fans_sigma(a, v, terms):
    sig=set()
    for p0 in all_simple(a, v, terms[0], set(terms)-{terms[0]}):
        u0=set(p0[1:])
        for p1 in all_simple(a, v, terms[1], (u0-{v}) | {terms[2]}):
            u1=u0|set(p1[1:])
            for p2 in all_simple(a, v, terms[2], (u1-{v})):
                tr=((len(p0)-1)%2,(len(p1)-1)%2,(len(p2)-1)%2)
                if len(set(tr))==1: sig.add(tr[0])
    return sig
def dist(a,s):
    d={s:0};q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w not in d: d[w]=d[u]+1;q.append(w)
    return d
def path_parities(a,v,t):
    return {(len(p)-1)%2 for p in all_simple(a,v,t,set())}
def on_odd_cycle(a,v):
    for w in a[v]:
        E2={frozenset((x,y)) for x in a for y in a[x]}-{frozenset((v,w))}
        a2=adj_of(E2)
        for p in all_simple(a2,v,w,set()):
            if (len(p)-1+1)%2==1: return True
    return False

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 191)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 14):
    a,terms=make_bridge(rng.choice([8,10]),rng)
    if a is None: continue
    D=[dist(a,t) for t in terms]
    st['bridges']+=1
    any0=False
    for v in a:
        if v in terms: continue
        S=fans_sigma(a,v,terms)
        h1=all(D[m][v]%2==0 for m in range(3))
        h2=all(len(path_parities(a,v,terms[m]))==2 for m in range(3))
        h3=on_odd_cycle(a,v)
        has0=(0 in S); has1=(1 in S)
        if has0: any0=True
        st['sigma0=%s  H1(dist even)=%s'%(has0,h1)]+=1
        st['sigma0=%s  H2(parity-flexible)=%s'%(has0,h2)]+=1
        st['sigma0=%s  H3(on odd cycle)=%s'%(has0,h3)]+=1
    st['bridge HAS an even balanced fan' if any0 else '*** NO even balanced fan ***']+=1
for k in sorted(st): print('  %-42s %s'%(k,st[k]))
