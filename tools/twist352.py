"""twist352.py -- large ODD GIRTH with tau >= 2: does the Terminal Cycle Theorem survive?

348.3 suggested large odd girth might force small tau, putting gap (i) inside the range 308
already closed.  That needs testing, not assuming: two DISJOINT long odd cycles give large
odd girth AND tau >= 2 at once.
Construction: random BIPARTITE cubic graph of girth >= 6; apply TWO distant TWISTS - replace
edges u1v1, u2v2 by u1u2 and v1v2, which joins same-side vertices and creates odd cycles.
Two twists far apart give tau = 2 with the odd girth set by their separation.
Then subdivide three far-apart edges and test the theorem.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def rand_bip_cubic(m,rng,tries=3000):
    for _ in range(tries):
        X=list(range(m)); Y=list(range(m,2*m))
        perms=[rng.sample(Y,m) for _ in range(3)]
        E=set(); ok=True
        for k in range(3):
            for i,x in enumerate(X):
                e=frozenset((x,perms[k][i]))
                if e in E: ok=False; break
                E.add(e)
            if not ok: break
        if not ok: continue
        a={v:set() for v in X+Y}
        for e in E:
            u,v=tuple(e); a[u].add(v); a[v].add(u)
        if any(len(a[v])!=3 for v in a): continue
        if girth(a)>=6: return a,E
    return None,None

def odd_girth(a):
    best=10**9
    for s in a:
        d={(s,0):0}; q=deque([(s,0)])
        while q:
            u,p=q.popleft()
            for w in a[u]:
                st=(w,1-p)
                if st not in d: d[st]=d[(u,p)]+1; q.append(st)
        if (s,1) in d: best=min(best,d[(s,1)])
    return best
def tau_edge(E,cap=3):
    El=sorted(E,key=lambda s:sorted(map(str,s)))
    for k in range(cap+1):
        for S in itertools.combinations(El,k):
            aa=adjof(set(El)-set(S))
            if bipartite(aa): return k
    return cap+1
def adjof(E):
    a={}
    for e in E:
        u,v=tuple(e); a.setdefault(u,set()).add(v); a.setdefault(v,set()).add(u)
    return a
def bipartite(a):
    col={}
    for s in a:
        if s in col: continue
        col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w not in col: col[w]=1-col[u]; q.append(w)
                elif col[w]==col[u]: return False
    return True
def both_parities(a,drop,s,t,rng,K=300):
    if s in drop or t in drop: return None
    par=set()
    for _ in range(K):
        p=[s]; on={s}|set(drop)
        while p[-1]!=t:
            nb=[w for w in a[p[-1]] if w not in on]
            if not nb: break
            w=rng.choice(nb); p.append(w); on.add(w)
        if p[-1]==t:
            par.add((len(p)-1)%2)
            if len(par)==2: return True
    return False if par else None

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 4201)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 20):
    a0,E0=rand_bip_cubic(16,rng)
    if a0 is None: st['bipartite gen failed']+=1; continue
    El=[tuple(e) for e in E0]; D0={v:dist(a0,v) for v in a0}
    done=None
    for _ in range(400):
        (u1,v1),(u2,v2)=rng.sample(El,2)
        if u1>=16: u1,v1=v1,u1
        if u2>=16: u2,v2=v2,u2
        if min(D0[x][y] for x in (u1,v1) for y in (u2,v2))<4: continue
        E=set(E0)-{frozenset((u1,v1)),frozenset((u2,v2))}
        E|={frozenset((u1,u2)),frozenset((v1,v2))}
        a=adjof(E)
        if any(len(a[v])!=3 for v in a): continue
        if girth(a)<5: continue
        if bipartite(a): continue
        done=(a,E); break
    if done is None: st['twist failed']+=1; continue
    a,E=done
    og=odd_girth(a); t=tau_edge(E)
    st['odd girth %2d, tau %d'%(og,t)]+=1
    if t<2: continue
    # subdivide three far-apart edges
    El2=[tuple(e) for e in E]; DD={v:dist(a,v) for v in a}
    br=None
    for _ in range(500):
        tri=rng.sample(El2,3); ok=True
        for (x1,y1),(x2,y2) in itertools.combinations(tri,2):
            if min(DD[p][q] for p in (x1,y1) for q in (x2,y2))<3: ok=False; break
        if not ok: continue
        E2=set(E); terms=[]
        for k,(u,v) in enumerate(tri):
            tt='b%d'%k; terms.append(tt)
            E2.discard(frozenset((u,v))); E2|={frozenset((u,tt)),frozenset((tt,v))}
        ab=adjof(E2); Db={x:dist(ab,x) for x in terms}
        if min(Db[terms[i]][terms[j]] for i,j in [(0,1),(0,2),(1,2)])>=5: br=(ab,terms); break
    if br is None: st['  bridge extraction failed']+=1; continue
    ab,terms=br
    seen=set(); applies=False
    for _ in range(6000):
        c=rand_cycle(ab,rng)
        if not c or len(c)%2 or not any(t2 in c for t2 in terms): continue
        k=frozenset(c)
        if k in seen: continue
        seen.add(k)
        Ds=set(c); n=len(c); offs=[t2 for t2 in terms if t2 not in Ds]
        for k2 in range(n):
            x,xp=c[k2],c[(k2+1)%n]
            ys=[y for y in ab[x] if y not in Ds]; yps=[y for y in ab[xp] if y not in Ds]
            if not(ys and yps): continue
            for bi in offs:
                if both_parities(ab,Ds,bi,ys[0],rng) and both_parities(ab,Ds,bi,yps[0],rng):
                    applies=True; break
            if applies: break
        if applies or len(seen)>60: break
    st['  tau=%d og=%d -> THEOREM APPLIES: %s'%(t,og,applies)]+=1
for k in sorted(st): print('  %s  %s'%(k.ljust(46),st[k]))
