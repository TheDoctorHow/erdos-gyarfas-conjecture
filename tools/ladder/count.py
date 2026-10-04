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

def escapes(adj,Zset,start,terms):
    seen={start}; q=deque([start])
    while q:
        u=q.popleft()
        for w in adj[u]:
            if w in Zset or w in seen: continue
            seen.add(w); q.append(w)
            if w in terms: return True
    return False

def survey(adj,b1,b2,terms,cap=600000):
    """count cycles through b1,b2 by (arc-parity class), and how many escape"""
    x,y=adj[b1]; cnt=0
    tot=Counter(); esc=0; mixed=0
    stack=[(x,1<<x|1<<b1,(x,))]
    while stack:
        u,mask,pth=stack.pop()
        cnt+=1
        if cnt>cap: return None
        for w in adj[u]:
            if w==y:
                cyc=(b1,)+pth+(y,)
                if b2 not in cyc: continue
                L=len(cyc); a=cyc.index(b2)
                pa,pb=a%2,(L-a)%2
                key="even-even" if (pa==0 and pb==0) else ("odd-odd" if (pa==1 and pb==1) else "mixed")
                tot[key]+=1
                if key=="mixed": mixed+=1
                if key=="even-even":
                    Z=set(cyc)
                    if any(t not in Z for t in terms):
                        for st in (x,y):
                            if st in Z and escapes(adj,Z,st,terms): esc+=1; break
            elif not (mask>>w & 1):
                stack.append((w,mask|1<<w,pth+(w,)))
    return tot,esc,mixed

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 101)
    ee=[]; eepar=Counter(); escpar=Counter(); nmixed=0; npairs=0; zero_ee=0
    for nH in (18,20):
        for trial in range(60):
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
                adj,term=build_B(nH,E,M); terms=set(term)
                for b1,b2 in itertools.combinations(term,2):
                    r=survey(adj,b1,b2,terms)
                    if r is None: continue
                    tot,esc,mixed=r
                    npairs+=1
                    ee.append(tot["even-even"])
                    eepar[tot["even-even"]%2]+=1
                    escpar[esc%2]+=1
                    if tot["even-even"]==0: zero_ee+=1
                    if mixed>0: nmixed+=1
    ee.sort()
    print("terminal pairs surveyed:",npairs)
    print("  pairs with ZERO both-arcs-even cycles :",zero_ee)
    print("  pairs with a MIXED-parity cycle       :",nmixed)
    if ee:
        print("  #both-arcs-even cycles per pair: min=%d median=%d max=%d"%(ee[0],ee[len(ee)//2],ee[-1]))
    print("\nparity of the COUNT of both-arcs-even cycles:")
    for k in sorted(eepar): print("   count = %d mod 2 -> %d pairs"%(k,eepar[k]))
    print("parity of the COUNT of those that also escape:")
    for k in sorted(escpar): print("   count = %d mod 2 -> %d pairs"%(k,escpar[k]))

main()
