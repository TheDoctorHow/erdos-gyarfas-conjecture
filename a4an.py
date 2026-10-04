import sys,collections
def wins(spec,c):
    return [a for a in range(3,c//2+1) if all(L in spec for L in range(a,2*a+1) if L%2==0)]
rows=[]
for fn in sys.argv[1:]:
    for line in open(fn):
        p=line.split()
        if len(p)<5: continue
        g6=p[0]; n=int(p[1][2:]); bip=int(p[2][4:]); spec=int(p[3][5:],16); unk=int(p[4][4:],16)
        assert unk==0, ("UNDECIDED",g6)
        S=set(L for L in range(3,n+1) if (spec>>L)&1)
        if not S: continue
        c=max(S); W=wins(S,c)
        rows.append((g6,n,bip,c,S,W))
print("total graphs:",len(rows))
# ---- (a) bipartite
bp=[r for r in rows if r[2]==1]
print("bipartite:",len(bp))
noW=[r for r in bp if not r[5]]
print("  bipartite with NO DWL window (a>=3, 2a<=c):",len(noW))
for r in sorted(noW,key=lambda r:-r[3])[:12]:
    print("   ",r[0],"n=",r[1],"circ=",r[3],"spec=",sorted(r[4]))
# ratio maxW/(c/2)
def ratio(r):
    return (max(r[5])/(r[3]/2)) if r[5] else 0.0
print("  most fragile bipartite by maxW/(c/2):")
for r in sorted(bp,key=lambda r:(ratio(r),-r[3]))[:10]:
    print("   ",r[0],"n=",r[1],"circ=",r[3],"maxW=",(max(r[5]) if r[5] else None),
          "ratio=%.3f"%ratio(r),"spec=",sorted(r[4]))
# ---- (c) sparsity of the EVEN part
print("\n== even-spectrum sparsity ==")
worst=[]
for r in rows:
    g6,n,bip,c,S,W=r
    ev=[L for L in range(4,c+1,2)]
    present=[L for L in ev if L in S]
    missing=[L for L in ev if L not in S]
    # longest run of consecutive MISSING even values
    run=best=0
    for L in ev:
        if L not in S: run+=1; best=max(best,run)
        else: run=0
    worst.append((best,len(missing),len(present),r))
worst.sort(key=lambda t:(-t[0],-t[1]))
print(" largest run of consecutive missing EVEN lengths inside [4,circ]:")
for t in worst[:12]:
    b,nm,np_,r=t
    print("   run=%d missing=%d present=%d  %s n=%d bip=%d circ=%d spec=%s W=%s"%(b,nm,np_,r[0],r[1],r[2],r[3],sorted(r[4]),r[5]))
# no-window graphs overall
nw=[r for r in rows if not r[5]]
print("\n graphs with EMPTY DWL window set:",len(nw))
byc=collections.Counter(r[3] for r in nw)
print("  their circumferences:",dict(sorted(byc.items())))
for r in sorted(nw,key=lambda r:-r[3])[:10]:
    print("   ",r[0],"n=",r[1],"bip=",r[2],"circ=",r[3],"spec=",sorted(r[4]))
# smallest maxW relative to c overall
print("\n most fragile overall by maxW/(c/2) among circ>=12:")
cand=[r for r in rows if r[3]>=12]
for r in sorted(cand,key=lambda r:(ratio(r),-r[3]))[:12]:
    print("   ",r[0],"n=",r[1],"bip=",r[2],"circ=",r[3],"maxW=",(max(r[5]) if r[5] else None),
          "ratio=%.3f"%ratio(r),"spec=",sorted(r[4]))
