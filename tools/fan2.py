"""fan2.py -- what actually witnesses the failure of rigidity in the (N3) graphs?

For each H (bipartite cubic minus two independent edges, four degree-2 vertices two per
class) and the pair (x,b_2) in DIFFERENT classes, rigidity says all x-b_2 paths agree
mod 4.  It fails.  This finds the SIMPLEST certificate:

  a pair of x-b_2 paths P,P' with |P| =/= |P'| mod 4 whose difference is a SINGLE CYCLE Z.

For such a pair, P meets Z in a single arc, P' in the complementary arc, and
|Z| =/= 2|P n Z| mod 4 -- exactly the parity test of 302.3 firing.  If a single-cycle
certificate always exists, (FAN-2) is the right target; the structure of the smallest Z
is then reported.
"""
import itertools, collections
from collections import deque
exec(open('tools/n3search.py').read().split('def heawood()')[0].split('"""',2)[2])

def heawood():
    E = set()
    for i in range(14): E.add(frozenset((i,(i+1)%14)))
    for i in range(0,14,2): E.add(frozenset((i,(i+5)%14)))
    return E

def eset(p): return frozenset(frozenset((p[i],p[i+1])) for i in range(len(p)-1))

def ncomp(D):
    a = adj(list(D)); seen=set(); c=0
    for s in a:
        if s in seen: continue
        c+=1; q=deque([s]); seen.add(s)
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in seen: seen.add(w); q.append(w)
    return c

st = collections.Counter(); zlens = collections.Counter(); arcinfo = collections.Counter()
E0 = heawood(); El = sorted(E0, key=lambda f: sorted(map(str,f)))
for e1, e2 in itertools.combinations(El, 2):
    if e1 & e2: continue
    H = set(E0)-{e1,e2}; aH = adj(H)
    if not two_conn(aH) or girth(aH) < 6: continue
    cH = bip(aH)
    if cH is None: continue
    d2 = [v for v in aH if len(aH[v])==2]
    if len(d2)!=4: continue
    X=[v for v in d2 if cH[v]==0]; Y=[v for v in d2 if cH[v]==1]
    if len(X)!=2 or len(Y)!=2: continue
    for x in X:
        for b2 in Y:
            P = [p for p in paths(aH,x,b2)]
            Es = [(eset(p), len(p)-1) for p in P]
            st['pairs'] += 1
            # does rigidity fail at all?
            if len({L%4 for _,L in Es}) == 1:
                st['rigid (unexpected)'] += 1; continue
            st['rigidity fails'] += 1
            best = None
            for i in range(len(Es)):
                for j in range(i+1,len(Es)):
                    if (Es[i][1]-Es[j][1]) % 4 == 0: continue
                    D = Es[i][0] ^ Es[j][0]
                    if ncomp(D) != 1: continue          # single cycle
                    if best is None or len(D) < len(best[0]):
                        arc = len(Es[i][0] & D)
                        best = (D, Es[i][1], Es[j][1], arc)
            if best is None:
                st['*** no single-cycle certificate ***'] += 1
            else:
                st['single-cycle certificate found'] += 1
                Z, L1, L2, arc = best
                zlens[len(Z)] += 1
                arcinfo['|Z|=%d  arc=%d  |Z|/2 mod2=%d  arc mod2=%d'
                        % (len(Z), arc, (len(Z)//2)%2, arc%2)] += 1
for k in sorted(st): print('  %-36s %s' % (k, st[k]))
print('  smallest certificate cycle length:', dict(sorted(zlens.items())))
for k in sorted(arcinfo): print('    %-44s %s' % (k, arcinfo[k]))
