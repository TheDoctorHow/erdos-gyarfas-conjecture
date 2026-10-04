"""extremal311.py -- the extremal family of 310.4, settled.

Equality in 310.3 forces K to be a path v_1..v_2m plus m CHERRIES, each adjacent to two
path vertices, the attachments forming a PERFECT MATCHING of the path, with every span
q-p == 2 (mod 4) and >= 6 (girth).  Spans are even, so m is even.

(1) For which m does such a matching exist?
(2) Verify the smallest one really is 2-connected, bipartite, girth 6, subcubic, has the
    predicted degree-2 set, and SATISFIES rigidity - i.e. the engine genuinely stalls.
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
def two_conn(a): return conn(a) and all(conn(a,v) for v in a)
def bip(a):
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
def paths(a,s,t):
    out=[];p=[s];on={s}
    def rec(u):
        for w in a[u]:
            if w==t: out.append(list(p)+[t])
            elif w not in on:
                on.add(w);p.append(w);rec(w);p.pop();on.discard(w)
    rec(s);return out

print('(1) does a valid cherry matching exist?   spans >= 6 and == 2 mod 4')
found={}
for m in range(2, 13, 2):
    V=list(range(1,2*m+1)); ok=None
    def match(rem, acc):
        global ok
        if ok: return
        if not rem: ok=list(acc); return
        p=rem[0]
        for q in rem[1:]:
            d=q-p
            if d>=6 and d%4==2:
                match([z for z in rem if z not in (p,q)], acc+[(p,q)])
                if ok: return
    ok=None; match(V,[])
    found[m]=ok
    print('   m=%2d (n=%2d, s=%2d): %s' % (m, 3*m, m+2, ok if ok else 'IMPOSSIBLE'))

print()
mm=min(k for k in found if found[k])
print('smallest extremal m = %d, so smallest extremal s = %d, n = %d'%(mm,mm+2,3*mm))
pairs=found[mm]
E=set()
for i in range(1,2*mm): E.add(frozenset(('v%d'%i,'v%d'%(i+1))))
for k,(p,q) in enumerate(pairs):
    E.add(frozenset(('w%d'%k,'v%d'%p))); E.add(frozenset(('w%d'%k,'v%d'%q)))
a=adj(E)
d2=sorted(v for v in a if len(a[v])==2)
print()
print('(2) the smallest extremal graph, pairs %s'%pairs)
print('   n=%d m=%d maxdeg=%d girth=%d 2-connected=%s bipartite=%s'
      %(len(a),len(E),max(len(a[v]) for v in a),girth(a),two_conn(a),bip(a) is not None))
print('   degree-2 vertices (%d): %s'%(len(d2),d2))
c=bip(a)
D=[v for v in d2 if v.startswith('w')]
dist=lambda s,t:(lambda d: d)(None)
def bfs(s):
    d={s:0};q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w not in d: d[w]=d[u]+1;q.append(w)
    return d
mind=min(bfs(u)[v] for u,v in itertools.combinations(D,2))
print('   min distance between non-end degree-2 vertices: %d  (need >= 3)'%mind)
L={len(p)-1 for p in paths(a,'v1','v%d'%(2*mm))}
print('   a-b path lengths: %s   mod 4: %s   RIGID: %s'
      %(sorted(L), sorted({l%4 for l in L}), len({l%4 for l in L})==1))
print('   |V| = %d, 3s-6 = %d  -> equality: %s'%(len(a),3*(mm+2)-6,len(a)==3*(mm+2)-6))
