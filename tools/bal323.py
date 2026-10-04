"""bal323.py -- structure of (BAL).

(BAL)  every bridge with tau >= 2 carries balanced fans of BOTH leg parities.
Stronger candidate:
(BAL+) some SINGLE vertex carries balanced fans of both parities.
If (BAL+) holds it is much the easier target: it asks that at one vertex all three legs can
be flipped in parity simultaneously, which is a statement about odd cycles reachable from
that vertex, not a global one.
Exhaustive fan enumeration; reports the population of each parity and the overlap.
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

def all_fans_sigma(a, v, terms):
    sig=set()
    for p0 in all_simple(a, v, terms[0], set(terms)-{terms[0]}):
        u0=set(p0[1:])
        for p1 in all_simple(a, v, terms[1], (u0-{v}) | {terms[2]}):
            u1=u0|set(p1[1:])
            for p2 in all_simple(a, v, terms[2], (u1-{v})):
                sig.add(((len(p0)-1)%2,(len(p1)-1)%2,(len(p2)-1)%2))
        if len(sig)==8: break
    return sig


rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 91)
st=collections.Counter(); rows=[]
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 30):
    a,terms=make_bridge(rng.choice([10,12,14]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    t=tau_free(E)
    if t<2: continue
    n0=set(); n1=set(); both=set(); nofan=0; alleight=0
    for v in a:
        if v in terms: continue
        S=all_fans_sigma(a,v,terms)
        if not S: nofan+=1; continue
        if len(S)==8: alleight+=1
        has0=(0,0,0) in S; has1=(1,1,1) in S
        if has0: n0.add(v)
        if has1: n1.add(v)
        if has0 and has1: both.add(v)
    st['bridges tau>=2']+=1
    st['(BAL)  holds' if (n0 and n1) else '(BAL)  FAILS']+=1
    st['(BAL+) holds' if both else '(BAL+) FAILS']+=1
    rows.append((len(a),t,len(n0),len(n1),len(both),nofan,alleight,len(a)-3))
print('   n  tau  #sigma0  #sigma1  #both  #no-fan  #all-8-triples  (of %s interior)'%'k')
for r in rows[:18]:
    print('  %2d   %d    %3d      %3d     %3d      %2d        %3d          %d'%r)
for k in sorted(st): print('  %-22s %s'%(k,st[k]))
tot=len(rows)
if tot:
    print('  mean fraction of interior vertices with a balanced fan of BOTH parities: %.2f'
          %(sum(r[4]/max(r[7],1) for r in rows)/tot))
    print('  mean fraction carrying ALL EIGHT parity triples: %.2f'
          %(sum(r[6]/max(r[7],1) for r in rows)/tot))
