"""const327.py -- applying the hypothesis-removal lens to 318.2 itself.

Sum the three (SC) relations of ONE balanced fan:
   (mu_1+mu_2) + (mu_1+mu_3) + (mu_2+mu_3) = 2(mu_1+mu_2+mu_3) == c_12 + c_13 + c_23 (mod 4)
and mu_1+mu_2+mu_3 == 3 sigma == sigma (mod 2), so 2(mu_1+mu_2+mu_3) == 2 sigma (mod 4):

   S := c_12 + c_13 + c_23 == 2 sigma  (mod 4),   i.e.  sigma == S/2  (mod 2).

So ONE balanced fan already pins sigma to an invariant of B.  The Balanced Fan Theorem is
the statement that a determined quantity cannot take two values.
Verified on bridges that actually SATISFY (SC) - i.e. where the constants c_ij exist.
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
def sc_consts(a,terms):
    """returns {pair: c_ij} if (SC) holds, else None"""
    C={}
    for i,j in itertools.combinations(range(3),2):
        L=[len(p)-1 for p in all_simple(a,terms[i],terms[j],set())]
        L=[x for x in L if x%2==0]
        if not L: return None
        cl={x%4 for x in L}
        if len(cl)>1: return None
        C[(i,j)]=cl.pop()
    return C
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

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 173)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 300):
    a,terms=make_bridge(rng.choice([8,10,12]),rng)
    if a is None: continue
    C=sc_consts(a,terms)
    if C is None:
        st['bridges where (SC) fails']+=1; continue
    st['bridges SATISFYING (SC)']+=1
    S=sum(C.values())%4
    st['  S = c12+c13+c23 mod 4 = %d'%S]+=1
    pred=(S//2)%2 if S%2==0 else None
    bad=0; nbal=0
    for v in a:
        if v in terms: continue
        for s in fans_sigma(a,v,terms):
            nbal+=1
            if pred is None or s!=pred: bad+=1
    st['  balanced fans matching sigma=S/2' if bad==0
       else '  *** MISMATCH ***']+=1
    st['  bridges with no balanced fan at all']+=1 if nbal==0 else 0
for k in sorted(st): print('  %-44s %s'%(k,st[k]))
