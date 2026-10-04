import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools, sys
from collections import deque, Counter
base=_HERE+""
exec(open(base+"surv.py").read().split("def main")[0])

def build_B(nH,E,M):
    Ms=set(M); adj=[[] for _ in range(nH+len(M))]; term=[]; nxt=nH
    for (a,b) in E:
        if (a,b) in Ms:
            adj[a].append(nxt); adj[nxt].append(a)
            adj[b].append(nxt); adj[nxt].append(b); term.append(nxt); nxt+=1
        else:
            adj[a].append(b); adj[b].append(a)
    return adj,term

def oddcycles(adj,N,cap=400000):
    out=set(); c=0
    for s in range(N):
        stack=[(s,(s,),1<<s)]
        while stack:
            u,p,mask=stack.pop()
            c+=1
            if c>cap: return list(out)
            for w in adj[u]:
                if w==s and len(p)>=3 and len(p)%2==1: out.add(frozenset(p))
                elif w>s and not (mask>>w & 1): stack.append((w,p+(w,),mask|1<<w))
    return list(out)

def nu_odd(Z):
    best=0
    def rec(i,used,cnt):
        nonlocal best
        best=max(best,cnt)
        if best>=3 or i>=len(Z): return
        for j in range(i,len(Z)):
            if not (Z[j]&used): rec(j+1,used|Z[j],cnt+1)
    rec(0,frozenset(),0); return best

def two_disj_even(adj,N,s,t,cap=120000):
    ps=[]; stack=[(s,1<<s,(s,))]; c=0
    while stack:
        u,mask,p=stack.pop()
        c+=1
        if c>cap: break
        for w in adj[u]:
            if w==t:
                if len(p)%2==0: ps.append(p)
            elif not (mask>>w & 1): stack.append((w,mask|1<<w,p+(w,)))
    for i in range(len(ps)):
        for j in range(i+1,len(ps)):
            if not (set(ps[i][1:]) & set(ps[j][1:-1])): return True
    return False

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 202)
    joint=Counter()
    for nH in (18,20):
        for trial in range(70):
            g=randcubic(nH,rng)
            if not g: continue
            E,adj0=g
            D=[dist(adj0,nH,v) for v in range(nH)]
            def edist(e,f): return min(D[a][b] for a in e for b in f)
            cand=[t for t in itertools.combinations(E,3)
                  if all(edist(t[i],t[j])>=3 for i in range(3) for j in range(i+1,3))]
            if not cand: continue
            rng.shuffle(cand)
            for M in cand[:2]:
                adj,term=build_B(nH,E,M); N=nH+len(M)
                nu=nu_odd(oddcycles(adj,N))
                for b1,b2 in itertools.combinations(term,2):
                    ok=two_disj_even(adj,N,b1,b2)
                    joint[(min(nu,3),ok)]+=1
    print("(nu_odd(B) capped at 3, pair has two disjoint EVEN paths?) -> count")
    for k,v in sorted(joint.items()): print("   nu=%s  L2**=%s -> %d"%(k[0],k[1],v))

main()
