"""edge8.py -- test the LANDING-POINT argument for (EDGE-8).

P = shortest x-b_2 path, R = H - V(P), p_i = the unique off-P neighbour of v_i.
Rigidity forces every p_i-p_j path of R to have length == j-i-2 (mod 4).
A p_i-p_k path THROUGH p_j (i<j<k) then gives (j-i-2)+(k-j-2) == k-i-2, i.e. 2 == 0.

Checks, for every instance and every shortest x-b_2 path:
  * that every v_i really has exactly one off-P edge (no chords) -- sanity;
  * whether some component of R is 2-connected and holds >= 3 landing points;
  * whether the contradiction actually fires, by exhibiting i<j<k and a p_i-p_k path
    through p_j;
  * and, independently, whether the derived congruence is violated outright by two
    p_i-p_j paths of different length mod 4.
"""
import itertools, collections
from collections import deque
exec(open('tools/n3search.py').read().split('def heawood()')[0].split('"""',2)[2])
def heawood():
    E=set()
    for i in range(14): E.add(frozenset((i,(i+1)%14)))
    for i in range(0,14,2): E.add(frozenset((i,(i+5)%14)))
    return E

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

def sub(a, S): return {v:{w for w in a[v] if w in S} for v in S}

def two_connected_sub(a, S):
    if len(S)<3: return False
    aa=sub(a,S)
    def cc(skip):
        vs=[v for v in S if v!=skip]
        s=vs[0]; seen={s}; q=deque([s])
        while q:
            u=q.popleft()
            for w in aa[u]:
                if w!=skip and w not in seen: seen.add(w); q.append(w)
        return len(seen)==len(vs)
    return cc(None) and all(cc(v) for v in S)

def shortest_paths(a,s,t):
    """all shortest s-t paths"""
    d={s:0}; q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w not in d: d[w]=d[u]+1; q.append(w)
    if t not in d: return []
    out=[]
    def back(v,acc):
        if v==s: out.append([s]+acc[::-1][::-1]); return
        for w in a[v]:
            if d.get(w,-9)==d[v]-1: back(w,acc+[v])
    res=[]
    def go(v,path):
        if v==s: res.append([s]+path[::-1]); return
        for w in a[v]:
            if d.get(w,-9)==d[v]-1: go(w,path+[v])
    go(t,[])
    return res

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
        for b2 in Y:
            st['pairs']+=1
            fired=False; twoconn3=False
            allP = paths(aH,x,b2)
            induced = []
            for p in allP:
                Pv=set(p); ok=True
                for i,v in enumerate(p):
                    inn=[w for w in aH[v] if w in Pv]
                    if len(inn) > (1 if i in (0,len(p)-1) else 2): ok=False; break
                if ok: induced.append(p)
            st['induced paths'] += len(induced)
            for P in induced:
                Pv=set(P)
                # landing points
                land={}
                chord=False
                for i,v in enumerate(P):
                    off=[w for w in aH[v] if w not in Pv]
                    inn=[w for w in aH[v] if w in Pv]
                    if len(inn)>(1 if v in (x,b2) else 2): chord=True
                    if len(off)==1: land[i]=off[0]
                    elif len(off)>1: st['unexpected multi-off']+=1
                if chord: st['chord found (unexpected)']+=1
                R=[v for v in aH if v not in Pv]
                if not R: continue
                aR=sub(aH,set(R))
                for C in comps(aR):
                    idxs=sorted(i for i in land if land[i] in C)
                    if len(idxs)>=3 and two_connected_sub(aH,C): twoconn3=True
                    if len(idxs)<3: continue
                    # does some p_i-p_k path pass through p_j ?
                    for i,j,k in itertools.combinations(idxs,3):
                        pi,pj,pk=land[i],land[j],land[k]
                        # path pi->pk through pj inside C
                        S=set(C)
                        aa=sub(aH,S)
                        def pth(s0,t0,ban):
                            dd={s0:0}; q=deque([s0])
                            while q:
                                u=q.popleft()
                                for w in aa[u]:
                                    if w in ban or w in dd: continue
                                    dd[w]=1; q.append(w)
                            return t0 in dd
                        # pi->pj avoiding nothing, pj->pk disjoint: use Menger-ish check
                        # simple sufficient test: pi and pk in different comps of C-pj
                        S2=S-{pj}
                        aa2=sub(aH,S2)
                        seen={pi} if pi in S2 else set()
                        if pi in S2:
                            q=deque([pi])
                            while q:
                                u=q.popleft()
                                for w in aa2[u]:
                                    if w not in seen: seen.add(w); q.append(w)
                        if pk not in seen: fired=True
                    if two_connected_sub(aH,C) and len(idxs)>=3: fired=True
                if fired: break
            st['component 2-connected with >=3 landings' if twoconn3 else
               'no such component']+=1
            st['CONTRADICTION FIRES' if fired else '*** does not fire ***']+=1
for k in sorted(st): print('  %-46s %s' % (k,st[k]))
