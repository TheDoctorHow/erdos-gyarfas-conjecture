"""both306.py -- (SC) gives BOTH rigidities: (x,b_2) and (y,b_3) (299.2a).
Contradicting EITHER kills the configuration.  Does the three-term collapse of 305.2 fire
for at least one of the two pairs, in every instance?

In the (N3) configuration X = {x, b_3} and Y = {y, b_2}, so once x and b_2 are chosen the
second rigid pair is (Y minus b_2, X minus x).
"""
import itertools, collections
from collections import deque
exec(open('tools/n3search.py').read().split('def heawood()')[0].split('"""',2)[2])
def heawood():
    E=set()
    for i in range(14): E.add(frozenset((i,(i+1)%14)))
    for i in range(0,14,2): E.add(frozenset((i,(i+5)%14)))
    return E
def sub(a,S): return {v:{w for w in a[v] if w in S} for v in S}
def comps(a):
    seen=set(); out=[]
    for s in a:
        if s in seen: continue
        blk={s}; q=deque([s]); seen.add(s)
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in seen: seen.add(w); blk.add(w); q.append(w)
        out.append(blk)
    return out
def two_connected_sub(a,S):
    if len(S)<3: return False
    aa=sub(a,S)
    def cc(skip):
        vs=[v for v in S if v!=skip]
        if not vs: return True
        s=vs[0]; seen={s}; q=deque([s])
        while q:
            u=q.popleft()
            for w in aa[u]:
                if w!=skip and w not in seen: seen.add(w); q.append(w)
        return len(seen)==len(vs)
    return cc(None) and all(cc(v) for v in S)

def collapse_fires(aH, s, t):
    """does the three-term collapse of 305.2 fire for the pair (s,t)?"""
    for p in paths(aH, s, t):
        Pv=set(p); ok=True
        for i,v in enumerate(p):
            inn=[w for w in aH[v] if w in Pv]
            if len(inn)>(1 if i in (0,len(p)-1) else 2): ok=False; break
        if not ok: continue
        land={}
        for i,v in enumerate(p):
            off=[w for w in aH[v] if w not in Pv]
            if len(off)==1: land[i]=off[0]
        R=set(aH)-Pv
        if len(R)<3: continue
        for C in comps(sub(aH,R)):
            idxs=sorted(i for i in land if land[i] in C)
            if len(idxs)<3: continue
            if two_connected_sub(aH,C): return True
            # or: p_j separates p_i from p_k inside C
            for i,j,k in itertools.combinations(idxs,3):
                pj=land[j]
                S2=C-{pj}
                if not S2: continue
                aa2=sub(aH,S2)
                pi,pk=land[i],land[k]
                if pi not in S2 or pk not in S2: continue
                seen={pi}; q=deque([pi])
                while q:
                    u=q.popleft()
                    for w in aa2[u]:
                        if w not in seen: seen.add(w); q.append(w)
                if pk not in seen: return True
    return False

st=collections.Counter()
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
    for x in X:
        b3=[v for v in X if v!=x][0]
        for b2 in Y:
            y=[v for v in Y if v!=b2][0]
            st['instances']+=1
            f1=collapse_fires(aH,x,b2)
            f2=collapse_fires(aH,y,b3)
            st['pair1 (x,b2) fires']+=f1
            st['pair2 (y,b3) fires']+=f2
            if f1 or f2: st['AT LEAST ONE FIRES']+=1
            else: st['*** NEITHER FIRES ***']+=1
for k in sorted(st): print('  %-34s %s'%(k,st[k]))
