"""probe303.py -- which violations of the parity test are available?

Rigidity forces, for every cycle Z met by an x-b_2 path P in a single arc A:
        |A| == |Z|/2  (mod 2).
So a violation is any clean (P,Z) with |A| of the wrong parity.  Two families:
   (a) |Z| = 2 mod 4 (e.g. a 6-cycle) met in an EVEN arc      -- 303.2's certificate
   (b) |Z| = 0 mod 4 (e.g. an 8-cycle) met in an ODD arc
and the sharpest case of (b) is arc = 1: a single EDGE of P lying on a clean 8-cycle,
which by 303.4 rigidity forbids outright.
"""
import itertools, collections
from collections import deque
exec(open('tools/n3search.py').read().split('def heawood()')[0].split('"""',2)[2])
def heawood():
    E=set()
    for i in range(14): E.add(frozenset((i,(i+1)%14)))
    for i in range(0,14,2): E.add(frozenset((i,(i+5)%14)))
    return E
def eset(p): return frozenset(frozenset((p[i],p[i+1])) for i in range(len(p)-1))

def all_cycles(a, maxlen=14):
    out=set(); vs=sorted(a,key=str); idx={v:i for i,v in enumerate(vs)}
    for i,s in enumerate(vs):
        p=[s]; on={s}
        def rec(u):
            if len(p)>maxlen: return
            for w in a[u]:
                if idx[w]<i: continue
                if w==s and len(p)>=3: out.add(frozenset(p))
                elif w not in on:
                    on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
        rec(s)
    return out

st=collections.Counter(); minarc=collections.Counter()
E0=heawood(); El=sorted(E0,key=lambda f:sorted(map(str,f)))
for e1,e2 in itertools.combinations(El,2):
    if e1&e2: continue
    H=set(E0)-{e1,e2}; aH=adj(H)
    if not two_conn(aH) or girth(aH)<6: continue
    cH=bip(aH)
    if cH is None: continue
    d2=[v for v in aH if len(aH[v])==2]
    if len(d2)!=4: continue
    X=[v for v in d2 if cH[v]==0]; Y=[v for v in d2 if cH[v]==1]
    if len(X)!=2 or len(Y)!=2: continue
    CY=[]
    for Zv in all_cycles(aH):
        Zed={frozenset((u,w)) for u in Zv for w in aH[u] if w in Zv}
        if len(Zed)==len(Zv): CY.append((Zv,Zed))
    for x in X:
        for b2 in Y:
            st['pairs']+=1
            P=[p for p in paths(aH,x,b2)]
            fam_a=fam_b=False; best=99
            for p in P:
                pe=eset(p); pv=set(p)
                ai_all={v:i for i,v in enumerate(p)}
                for Zv,Zed in CY:
                    common=[v for v in p if v in Zv]
                    if not common: continue
                    ed=pe&Zed
                    if not ed: continue
                    if len(common)!=len(ed)+1: continue
                    ii=[ai_all[v] for v in common]
                    if max(ii)-min(ii)!=len(ii)-1: continue
                    A=len(ed); Zl=len(Zed)
                    if A%2 != (Zl//2)%2:
                        best=min(best,A)
                        if Zl%4==2: fam_a=True
                        else: fam_b=True
            st['family (a) 2-mod-4 cycle, even arc' if fam_a else 'no family (a)']+=1
            st['family (b) 0-mod-4 cycle, odd arc' if fam_b else 'no family (b)']+=1
            minarc['smallest violating arc = %d'%best]+=1
for k in sorted(st): print('  %-44s %s' % (k,st[k]))
for k in sorted(minarc): print('   %-36s %s' % (k,minarc[k]))
