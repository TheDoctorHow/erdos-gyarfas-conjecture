import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools
base=_HERE+""
exec(open(base+"bicycle.py").read().split("def main")[0])

def allpaths(adj,s,t):
    out=[]; stack=[(s,{s},[])]
    while stack:
        u,seen,es=stack.pop()
        for w in adj[u]:
            e=(min(u,w),max(u,w))
            if w==t: out.append(es+[e])
            elif w not in seen: stack.append((w,seen|{w},es+[e]))
    return out

def main():
    rng=random.Random(2026)
    ok=0; bad=0; tested=0
    worst=None
    for n in range(8,15):
        for s in range(3,min(n,7)+1):
            for _ in range(400):
                g=gen(n,s,rng)
                if not g: continue
                E,adj=g
                if n-s<1: continue
                idx={e:i for i,e in enumerate(E)}
                mu=len(E)-n+1
                # bipartite?
                col={0:0}; q=[0]; bip=True
                while q:
                    u=q.pop()
                    for w in adj[u]:
                        if w not in col: col[w]=1-col[u]; q.append(w)
                        elif col[w]==col[u]: bip=False
                if bip: continue
                tested+=1
                u,v=0,1                      # two terminals
                P=allpaths(adj,u,v)
                if not P: continue
                ev=[p for p in P if len(p)%2==0]
                if not ev: continue
                P0=ev[0]; b0=0
                for e in P0: b0^=1<<idx[e]
                diffs=[]
                for p in ev[1:]:
                    b=0
                    for e in p: b^=1<<idx[e]
                    diffs.append(b^b0)
                dim=len(gauss(diffs))
                target=mu-1              # dim of even-weight part of cycle space
                if dim==target: ok+=1
                else:
                    bad+=1
                    if worst is None or target-dim>worst[0]:
                        worst=(target-dim,n,s,mu,dim,target)
    print("non-bipartite bridges tested          :",tested)
    print("span of even path-differences = C_0   :",ok)
    print("strictly smaller                      :",bad)
    if worst: print("  largest shortfall: n=%d s=%d mu=%d  dim=%d vs target %d"%worst[1:])

main()
