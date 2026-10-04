"""adv354.py -- the Terminal Cycle Theorem against the programme's OWN adversarial graph.

1.2 constructs a cubic graph of girth 5 whose only 5-cycle is an added pentagon and which has
NO 6-, 7- or 8-cycle: tools/pentagon70.g6, n = 70.  Subdividing three edges FAR FROM THE
PENTAGON makes (FIVE) fail outright - no terminal 6-cycle exists, and the shortest odd cycle
through those edges is at least 9, so the shortest even terminal cycle is at least 10.
Does the Terminal Cycle Theorem still apply?
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def read_g6(path):
    s=open(path).read().split()[0]
    if s.startswith('>>graph6<<'): s=s[10:]
    b=[ord(c)-63 for c in s]
    if b[0]==63:                       # '~' prefix: n encoded in next 3 bytes
        n=(b[1]<<12)|(b[2]<<6)|b[3]; rest=b[4:]
    else:
        n=b[0]; rest=b[1:]
    bits=''.join(format(x,'06b') for x in rest)
    a={v:set() for v in range(n)}
    k=0
    for j in range(1,n):
        for i in range(j):
            if k<len(bits) and bits[k]=='1': a[i].add(j); a[j].add(i)
            k+=1
    return a

def both_parities(a,drop,s,t,rng,K=400):
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

def shortest_even_terminal_cycle(a,terms):
    best=None
    for b in terms:
        nb=sorted(a[b],key=str)
        if len(nb)!=2: continue
        p,q=nb
        d={(p,0):0}; Q=deque([(p,0)])
        while Q:
            u,par=Q.popleft()
            for w in a[u]:
                if w==b: continue
                s=(w,1-par)
                if s not in d: d[s]=d[(u,par)]+1; Q.append(s)
        if (q,0) in d:
            L=d[(q,0)]+2
            if best is None or L<best[0]: best=(L,b)
    return best

G=read_g6('tools/pentagon70.g6')
print('  loaded pentagon70: n=%d, cubic=%s, girth=%d'
      %(len(G), all(len(G[v])==3 for v in G), girth(G)))
rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 5101)
# locate the unique pentagon
pent=None
for s in G:
    for p1 in G[s]:
        for p2 in G[p1]:
            if p2==s: continue
            for p3 in G[p2]:
                if p3 in (s,p1): continue
                for p4 in G[p3]:
                    if p4 in (s,p1,p2): continue
                    if s in G[p4]: pent=set((s,p1,p2,p3,p4))
                if pent: break
            if pent: break
        if pent: break
    if pent: break
print('  unique 5-cycle at %s'%sorted(pent))
DP={v:dist(G,v) for v in G}
El=[tuple(e) for e in {frozenset((u,w)) for u in G for w in G[u]}]
far=[e for e in El if min(DP[x][p] for x in e for p in pent)>=3]
print('  edges at distance >=3 from the pentagon: %d'%len(far))
ok=0; tried=0
for trial in range(40):
    tri=rng.sample(far,3); good=True
    for (u1,v1),(u2,v2) in itertools.combinations(tri,2):
        if min(DP[x][y] for x in (u1,v1) for y in (u2,v2))<3: good=False; break
    if not good: continue
    E2={frozenset((u,w)) for u in G for w in G[u]}; terms=[]
    for k,(u,v) in enumerate(tri):
        t='b%d'%k; terms.append(t)
        E2.discard(frozenset((u,v))); E2|={frozenset((u,t)),frozenset((t,v))}
    a={}
    for e in E2:
        u,v=tuple(e); a.setdefault(u,set()).add(v); a.setdefault(v,set()).add(u)
    Db={t:dist(a,t) for t in terms}
    if min(Db[terms[i]][terms[j]] for i,j in [(0,1),(0,2),(1,2)])<5: continue
    tried+=1
    r=shortest_even_terminal_cycle(a,terms)
    if r is None: print('  no even terminal cycle at all'); continue
    L,b=r
    seen=set(); applies=False
    for _ in range(20000):
        c=rand_cycle(a,rng)
        if not c or len(c)!=L or b not in c: continue
        k=frozenset(c)
        if k in seen: continue
        seen.add(k)
        Ds=set(c); n=len(c); offs=[t for t in terms if t not in Ds]
        for k2 in range(n):
            x,xp=c[k2],c[(k2+1)%n]
            ys=[y for y in a[x] if y not in Ds]; yps=[y for y in a[xp] if y not in Ds]
            if not(ys and yps): continue
            for bi in offs:
                if both_parities(a,Ds,bi,ys[0],rng) and both_parities(a,Ds,bi,yps[0],rng):
                    applies=True; break
            if applies: break
        if applies or len(seen)>25: break
    print('  bridge %d: shortest even terminal cycle = %2d  -> THEOREM APPLIES: %s'
          %(tried,L,applies))
    if applies: ok+=1
    if tried>=6: break
print('  RESULT: theorem applies in %d of %d adversarial bridges'%(ok,tried))
