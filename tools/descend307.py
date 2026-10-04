"""descend307.py -- does the block descent close the 84 resistant instances?

Handle C, landing points p_i,p_j.  Every cut vertex of C separates them, so C is a chain
of blocks; rigidity for (p_i,p_j) in C forces rigidity for each block's entry/exit pair.
Take a block K that is 2-connected (not a single edge) with entry/exit a,b and test
whether rigidity for (a,b) in K is VIOLATED outright: two a-b paths of different length
mod 4.
"""
import itertools, collections
from collections import deque
exec(open('tools/both306.py').read().split('st=collections.Counter()')[0].split('"""',2)[2])

def blocks_of(a):
    idx={};low={};par={};stk=[];out=[];C=[0]
    def dfs(u):
        idx[u]=low[u]=C[0];C[0]+=1
        for w in a[u]:
            if w not in idx:
                par[w]=u; stk.append((u,w)); dfs(w); low[u]=min(low[u],low[w])
                if low[w]>=idx[u]:
                    blk=[]
                    while True:
                        e=stk.pop(); blk.append(e)
                        if e==(u,w): break
                    out.append(blk)
            elif w!=par.get(u) and idx[w]<idx[u]:
                low[u]=min(low[u],idx[w]); stk.append((u,w))
    for s in a:
        if s not in idx: par[s]=None; dfs(s)
    return out

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
            if collapse_fires(aH,x,b2): continue
            st['resistant']+=1
            closed=False
            for p in paths(aH,x,b2):
                Pv=set(p); ok=True
                for i,v in enumerate(p):
                    inn=[w for w in aH[v] if w in Pv]
                    if len(inn)>(1 if i in (0,len(p)-1) else 2): ok=False;break
                if not ok: continue
                land={}
                for i,v in enumerate(p):
                    off=[w for w in aH[v] if w not in Pv]
                    if len(off)==1: land[i]=off[0]
                for C in comps(sub(aH,set(aH)-Pv)):
                    idxs=sorted(i for i in land if land[i] in C)
                    if len(idxs)!=2: continue
                    pi,pj=land[idxs[0]],land[idxs[1]]
                    if pi==pj: continue
                    aC=sub(aH,C)
                    st['handles examined']+=1
                    # cut vertices of C must separate p_i from p_j -- verify
                    for z in C:
                        if z in (pi,pj): continue
                        S=C-{z}; aa=sub(aH,S)
                        seen={pi}; q=deque([pi])
                        while q:
                            u=q.popleft()
                            for w in aa[u]:
                                if w not in seen: seen.add(w); q.append(w)
                        if len(seen)!=len(S) and pj in seen:
                            st['*** cut vertex NOT separating ***']+=1
                    # blocks of C, and rigidity inside each 2-connected block
                    for blk in blocks_of(aC):
                        vs={v for e in blk for v in e}
                        if len(vs)<3: continue
                        st['2-connected blocks']+=1
                        aK=sub(aH,vs)
                        # entry/exit of this block along p_i..p_j
                        ends=[v for v in vs if v in (pi,pj)] + \
                             [v for v in vs if any(w not in vs for w in aC[v])]
                        ends=list(dict.fromkeys(ends))
                        for u,v in itertools.combinations(ends,2):
                            L={len(q)-1 for q in paths(aK,u,v)}
                            if len({l%4 for l in L})>1:
                                closed=True; st['block rigidity VIOLATED']+=1
                                break
                        if closed: break
                    if closed: break
                if closed: break
            st['DESCENT CLOSES IT' if closed else '*** descent stalls ***']+=1
for k in sorted(st): print('  %-38s %s'%(k,st[k]))
