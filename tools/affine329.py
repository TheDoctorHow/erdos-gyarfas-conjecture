"""affine329.py -- is the achievable parity-triple set Sigma(v) an AFFINE SUBSPACE of F_2^3?

The parity flips of Lemma B (326.2) act on a fan's parity triple.  If the available flips
generated a group, Sigma(v) would be a COSET of a subgroup: its size a power of two, and
closed under x+y+z.  That would explain 321.2(iv) - few vertices realising all eight - and
turn (BAL) into a statement about which cosets occur.
Tested exhaustively: the size distribution of Sigma(v), and the coset property.
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
def sigma_set(a, v, terms):
    S=set()
    for p0 in all_simple(a, v, terms[0], set(terms)-{terms[0]}):
        u0=set(p0[1:])
        for p1 in all_simple(a, v, terms[1], (u0-{v}) | {terms[2]}):
            u1=u0|set(p1[1:])
            for p2 in all_simple(a, v, terms[2], (u1-{v})):
                S.add(((len(p0)-1)%2,(len(p1)-1)%2,(len(p2)-1)%2))
        if len(S)==8: break
    return S
def is_coset(S):
    if not S: return True
    for x in S:
        for y in S:
            for z in S:
                w=tuple((x[i]^y[i]^z[i]) for i in range(3))
                if w not in S: return False
    return True

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 211)
st=collections.Counter(); bad=[]
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 12):
    a,terms=make_bridge(rng.choice([8,10]),rng)
    if a is None: continue
    st['bridges']+=1
    for v in a:
        if v in terms: continue
        S=sigma_set(a,v,terms)
        st['|Sigma(v)| = %d'%len(S)]+=1
        ok=is_coset(S)
        st['Sigma is a coset' if ok else '*** NOT a coset ***']+=1
        if not ok and len(bad)<4: bad.append((len(a),len(S),sorted(S)))
for k in sorted(st): print('  %-30s %s'%(k,st[k]))
if bad:
    print('  non-coset examples (n, |Sigma|, Sigma):')
    for b in bad: print('    ',b)
