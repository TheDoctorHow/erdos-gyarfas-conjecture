"""unbal336.py -- the route 318.3 dismissed too early.

Two fans, at v and v', BOTH having mu_i == mu_j (mod 2) for the SAME pair (i,j).  Then
L_i+L_j and L'_i+L'_j are both EVEN connecting paths for that pair, so (SC) pins both:
     mu_i + mu_j == c_ij == mu'_i + mu'_j   (mod 4)
hence  e_i + e_j == 0 (mod 4),  where e_m = mu_m - mu'_m.
318.3 dismissed this because MOD 2 it is automatic.  But the relation is MOD 4.
If some pair of such fans has e_i + e_j == 2 (mod 4), (SC) FAILS - with NO parity
prescription on the third leg, and no balanced fan required.
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
def fans(a, v, terms):
    """all fans at v, as triples of leg LENGTHS"""
    out=[]
    for p0 in all_simple(a, v, terms[0], set(terms)-{terms[0]}):
        u0=set(p0[1:])
        for p1 in all_simple(a, v, terms[1], (u0-{v}) | {terms[2]}):
            u1=u0|set(p1[1:])
            for p2 in all_simple(a, v, terms[2], (u1-{v})):
                out.append((len(p0)-1,len(p1)-1,len(p2)-1))
    return out

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 307)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 6):
    a,terms=make_bridge(10,rng)
    if a is None: continue
    st['bridges']+=1
    F={v:fans(a,v,terms) for v in a if v not in terms}
    hit=False; witness=None
    verts=[v for v in F if F[v]]
    for i,j in [(0,1),(0,2),(1,2)]:
        for v1,v2 in itertools.combinations(verts,2):
            for m1 in F[v1]:
                if (m1[i]-m1[j])%2: continue
                for m2 in F[v2]:
                    if (m2[i]-m2[j])%2: continue
                    e=( (m1[i]-m2[i]) + (m1[j]-m2[j]) ) % 4
                    if e!=0:
                        hit=True; witness=(i,j,m1,m2,e); break
                if hit: break
            if hit: break
        if hit: break
    st['*** MOD-4 RELATION VIOLATED -> (SC) fails ***' if hit
       else 'relation holds everywhere (no contradiction)']+=1
    if witness: st['  example pair=(%d,%d) legs %s vs %s, e_i+e_j=%d mod 4'%witness]+=1
for k in sorted(st): print('  %-56s %s'%(k,st[k]))
