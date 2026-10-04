import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools, sys
from collections import deque
exec(open(_HERE+"surv.py").read().split("def main")[0])

def build_B(nH,E,M):
    """subdivide each matching edge once; return adj, terminal list"""
    Ms=set(M); adj=[[] for _ in range(nH+len(M))]; term=[]
    nxt=nH
    for (a,b) in E:
        if (a,b) in Ms:
            adj[a].append(nxt); adj[nxt].append(a)
            adj[b].append(nxt); adj[nxt].append(b)
            term.append(nxt); nxt+=1
        else:
            adj[a].append(b); adj[b].append(a)
    return adj,term

def violates_SC(adj,n,s,t,cap=400000):
    """return (True, l1, l2) if two even s-t paths differ mod 4"""
    seen_class=None; cnt=0
    stack=[(s,1<<s,0)]
    while stack:
        u,mask,d=stack.pop()
        cnt+=1
        if cnt>cap: return None
        for w in adj[u]:
            if w==t:
                L=d+1
                if L%2==0:
                    if seen_class is None: seen_class=L%4
                    elif L%4!=seen_class: return (L,seen_class)
            elif not (mask>>w & 1):
                stack.append((w,mask|1<<w,d+1))
    return False

def main():
    rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 4)
    checked=0; sc_holds=0; sc_fails=0; timeouts=0
    examples=[]
    for nH in (18,20,22):
        for trial in range(120):
            g=randcubic(nH,rng)
            if not g: continue
            E,adj=g
            D=[dist(adj,nH,v) for v in range(nH)]
            def edist(e,f): return min(D[a][b] for a in e for b in f)
            cand=[t for t in itertools.combinations(E,3)
                  if all(edist(t[i],t[j])>=3 for i in range(3) for j in range(i+1,3))]
            if not cand: continue
            rng.shuffle(cand)
            for M in cand[:4]:
                Ms=set(M)
                neg=[]
                for C in cycles_upto(adj,nH,14):
                    ce=cycedges(adj,C)
                    if len(ce)!=len(C): continue
                    if (len(C)+sum(1 for e in ce if e in Ms))%2==1: neg.append(frozenset(C))
                best=0
                def rec(i,used,cnt):
                    nonlocal best
                    if cnt>best: best=cnt
                    if best>=3 or i>=len(neg): return
                    for j in range(i,len(neg)):
                        if not (neg[j]&used): rec(j+1,used|neg[j],cnt+1)
                rec(0,frozenset(),0)
                if best>=3: continue          # already excluded by 223.3
                bad,term=build_B(nH,E,M)
                checked+=1
                bad_pair=None; to=False
                for x,y in itertools.combinations(term,2):
                    r=violates_SC(bad,nH+len(M),x,y)
                    if r is None: to=True; break
                    if r: bad_pair=(x,y,r); break
                if to: timeouts+=1
                elif bad_pair: sc_fails+=1
                else:
                    sc_holds+=1
                    if len(examples)<4: examples.append((nH,M,best))
    print("configs with nu_odd <= 2, tested for (SC):",checked)
    print("  (SC) FAILS anyway (some pair violates) :",sc_fails)
    print("  (SC) HOLDS  -> would be a counterexample:",sc_holds)
    print("  timed out                              :",timeouts)
    for e in examples: print("   HOLDS:",e)

main()
