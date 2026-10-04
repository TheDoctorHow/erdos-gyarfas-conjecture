"""ladder312.py -- can the STALL CASE of the descent ever be rigid?

Stall case: every handle of the chordless path is a tree, hence (all its vertices having
degree 2) a PATH.  So K = path v_0..v_L plus internally disjoint subdivided chords, each
path vertex in at most one chord, v_0 and v_L chord endpoints (they have degree 2).
Constraints: subcubic (automatic), bipartite (chord length == span mod 2), girth >= 6,
2-connected (every gap of the path covered by a chord).
Question: does ANY such K have all v_0-v_L paths congruent mod 4?
"""
import itertools, collections
from collections import deque

def adj(E):
    a={}
    for e in E:
        u,v=tuple(e); a.setdefault(u,set()).add(v); a.setdefault(v,set()).add(u)
    return a
def girth(a):
    best=10**9
    for s in a:
        d={s:0};par={s:None};q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in d: d[w]=d[u]+1;par[w]=u;q.append(w)
                elif w!=par[u]: best=min(best,d[u]+d[w]+1)
    return best
def conn(a,skip=None):
    vs=[v for v in a if v!=skip]
    if not vs: return True
    s=vs[0];seen={s};q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w!=skip and w not in seen: seen.add(w);q.append(w)
    return len(seen)==len(vs)
def two_conn(a): return len(a)>=3 and conn(a) and all(conn(a,v) for v in a)
def paths(a,s,t,cap=200000):
    out=[];p=[s];on={s}
    def rec(u):
        if len(out)>cap: return
        for w in a[u]:
            if w==t: out.append(len(p))
            elif w not in on:
                on.add(w);p.append(w);rec(w);p.pop();on.discard(w)
    rec(s);return out

def matchings(verts):
    """all partial matchings of verts (as ordered pairs p<q) covering the two ends"""
    verts=sorted(verts)
    res=[]
    def rec(rem, acc):
        if not rem: res.append(list(acc)); return
        # option: leave rem[0] unmatched (allowed unless it is an end)
        p=rem[0]
        if p not in (verts[0], verts[-1]):
            rec(rem[1:], acc)
        for i in range(1,len(rem)):
            q=rem[i]
            rec([z for z in rem[1:] if z!=q], acc+[(p,q)])
    rec(verts, [])
    return res

st=collections.Counter(); rigid=[]
for L in range(4, 11):
    V=list(range(L+1))
    for M in matchings(V):
        if not M: continue
        spans=[q-p for p,q in M]
        # chord lengths: same parity as span (bipartite), girth: l + span >= 6, l >= 1
        ranges=[]
        ok=True
        for (p,q) in M:
            sp=q-p
            opts=[l for l in range(1,5) if (l-sp)%2==0 and l+sp>=6]
            if not opts: ok=False;break
            ranges.append(opts)
        if not ok: continue
        for combo in itertools.product(*ranges):
            E=set()
            for i in range(L): E.add(frozenset(('v%d'%i,'v%d'%(i+1))))
            idc=0
            for (p,q),l in zip(M,combo):
                prev='v%d'%p
                for k in range(l-1):
                    nv='c%d_%d'%(idc,k); E.add(frozenset((prev,nv))); prev=nv
                E.add(frozenset((prev,'v%d'%q))); idc+=1
            a=adj(E)
            if max(len(a[v]) for v in a)>3: continue
            if girth(a)<6: continue
            if not two_conn(a): continue
            st['ladders tested']+=1
            Ls={l for l in paths(a,'v0','v%d'%L)}
            if len({x%4 for x in Ls})==1:
                st['*** RIGID LADDER ***']+=1
                if len(rigid)<4: rigid.append((L,M,combo,sorted(Ls)))
print('  ladders tested          ', st['ladders tested'])
print('  RIGID ones              ', st['*** RIGID LADDER ***'])
for r in rigid: print('   L=%d chords=%s lengths=%s  path lengths=%s'%r)
