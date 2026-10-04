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

def cycles_through(adj,N,b,cap=200000):
    """all cycles through degree-2 vertex b: b + (x-y path avoiding b)"""
    x,y=adj[b]; out=[]; cnt=0
    stack=[(x,1<<x|1<<b,(x,))]
    while stack:
        u,mask,pth=stack.pop()
        cnt+=1
        if cnt>cap: return None
        for w in adj[u]:
            if w==y: out.append(pth+(y,))
            elif not (mask>>w & 1): stack.append((w,mask|1<<w,pth+(w,)))
    return out

def escapes(adj,N,Zset,start,terms):
    """can we reach a terminal not in Z from `start`, avoiding Z except start?"""
    seen={start}; q=deque([start])
    while q:
        u=q.popleft()
        if u!=start and u in terms and u not in Zset: return True
        for w in adj[u]:
            if w in Zset or w in seen: continue
            seen.add(w); q.append(w)
            if w in terms: return True
    return False

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 21)
    tot=0; good=0; classes=Counter(); nogood=[]
    for nH in (18,20):
        for trial in range(120):
            g=randcubic(nH,rng)
            if not g: continue
            E,adj0=g
            D=[dist(adj0,nH,v) for v in range(nH)]
            def edist(e,f): return min(D[a][b] for a in e for b in f)
            cand=[t for t in itertools.combinations(E,3)
                  if all(edist(t[i],t[j])>=3 for i in range(3) for j in range(i+1,3))]
            if not cand: continue
            rng.shuffle(cand)
            for M in cand[:3]:
                adj,term=build_B(nH,E,M); N=nH+len(M); terms=set(term)
                tot+=1; hit=False
                for b in term:
                    cyc=cycles_through(adj,N,b)
                    if cyc is None: continue
                    x,y=adj[b]
                    for pth in cyc:
                        L=len(pth)+1          # |Z| = |path| edges +2 ; path has len(pth)-1 edges
                        Zset=set(pth)|{b}
                        # segment-1 realisation: start at x (or y), escape to another terminal
                        for st in (x,y):
                            if st not in Zset: continue
                            if escapes(adj,N,Zset,st,terms):
                                sigma=(L-2)%4
                                classes[(L%4,sigma)]+=1
                                if L%4==0: hit=True
                                break
                        if hit: break
                    if hit: break
                if hit: good+=1
                else: nogood.append(nH)
    print("surviving configurations tested                    :",tot)
    print("  with a terminal on a REALISABLE 0-mod-4 cycle    :",good)
    print("  without                                          :",len(nogood))
    print("\n(|Z| mod 4, shift mod 4) over realisable terminal cycles, top entries:")
    for k,v in sorted(classes.items()): print("   ",k,"->",v)

main()
