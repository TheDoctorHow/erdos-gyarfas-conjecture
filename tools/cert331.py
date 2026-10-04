"""cert331.py -- (c): can the certificate be contradicted, and by what?

If the target fails then for every v there is F_v with |F_v| <= 5 meeting every odd
v-T* path.  By Lemma B, B' - F_v must have a BIPARTITE v-T* block chain; if it stays
2-connected that means F_v BIPARTIZES B' with v and T* in the SAME colour class.
Define  beta(v) = least |F| with B'-F bipartite and v, T* the same colour.
The certificate then forces beta(v) <= 5 for EVERY v.  Measured against reality.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split(chr(34)*3,2)[2])

def bip_col(a):
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
def beta(E, v, cap=6):
    El=sorted(E,key=lambda s:sorted(map(str,s)))
    for k in range(cap+1):
        for S in itertools.combinations(El,k):
            c=bip_col(adj_of(set(El)-set(S)))
            if c is None: continue
            # same colour in EVERY component-consistent 2-colouring: flip components freely,
            # so v,T* can be made to agree unless they are in the same component with
            # different colours
            if c.get(v) is None or c.get('T*') is None: return k
            # components are independently flippable
            comp={}
            seen=set()
            for s in adj_of(set(El)-set(S)):
                if s in seen: continue
                st=[s]; seen.add(s)
                while st:
                    u=st.pop()
                    comp[u]=s
                    for w in adj_of(set(El)-set(S))[u]:
                        if w not in seen: seen.add(w); st.append(w)
            if comp.get(v)!=comp.get('T*') or c[v]==c['T*']: return k
    return cap+1

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 257)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 8):
    a,terms=make_bridge(8,rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    for x in terms: E.add(frozenset(('T*',x)))
    aB=adj_of(E)
    st['bridges']+=1
    mx=0
    for v in aB:
        if v=='T*': continue
        b=beta(E,v)
        mx=max(mx,b)
        st['beta(v) = %d'%b]+=1
    st['max beta over v = %d'%mx]+=1
    st['CERTIFICATE REFUTED (beta>5 somewhere)' if mx>5
       else 'certificate survives (beta<=5 everywhere)']+=1
for k in sorted(st): print('  %-44s %s'%(k,st[k]))
