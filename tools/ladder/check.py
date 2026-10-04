import os as _os; _HERE = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep
import random, itertools
base=_HERE+""
exec(open(base+"bicycle.py").read().split("def main")[0])

def oddcycles(adj,n):
    """all odd cycles as vertex frozensets"""
    out=set()
    for s in range(n):
        stack=[(s,[s],{s})]
        while stack:
            u,path,seen=stack.pop()
            for w in adj[u]:
                if w==s and len(path)>=3 and len(path)%2==1:
                    out.add(frozenset(path))
                elif w not in seen and w>s:
                    stack.append((w,path+[w],seen|{w}))
    return list(out)

def nu_odd(adj,n):
    Z=oddcycles(adj,n)
    best=0
    def rec(i,used,cnt):
        nonlocal best
        best=max(best,cnt)
        if i>=len(Z): return
        if cnt+(len(Z)-i)<=best: return
        for j in range(i,len(Z)):
            if not (Z[j]&used): rec(j+1,used|Z[j],cnt+1)
    rec(0,frozenset(),0)
    return best

def kerAD(adj,n):
    """dim ker(A+D) over F2, D=diag(deg mod 2)"""
    rows=[]
    for v in range(n):
        r=0
        for w in adj[v]: r^=1<<w
        if len(adj[v])%2: r^=1<<v
        rows.append(r)
    rk=len(gauss(rows))
    return n-rk

rng=random.Random(31); bad=0; tested=0; tight=0
for n in range(7,16):
    for s in range(3,min(n,8)+1):
        for _ in range(700):
            g=gen(n,s,rng)
            if not g: continue
            E,adj=g
            mu=len(E)-n+1
            cut,cyc=spaces(E,n); d=dim_inter(cut,cyc)
            # verify the matrix identity too
            if kerAD(adj,n)-1 != d:
                print("MATRIX IDENTITY FAILS", n,s); bad+=1; continue
            nu=nu_odd(adj,n)
            tested+=1
            if nu > mu-d:
                print("BOUND FAILS n=%d s=%d mu=%d d=%d nu=%d"%(n,s,mu,d,nu)); bad+=1
            elif nu == mu-d: tight+=1
print("graphs tested            :",tested)
print("violations of nu <= mu-d :",bad)
print("cases where bound is tight:",tight)
