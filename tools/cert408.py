"""cert408.py -- the certificate calculus of 408, audited.

Per bridge (tau>=2), per pair, per sampled even connecting path P with Z(P) the
component of b_k in B-P:
  (1) CRITERION: for every interior z of P, chain-criterion (odd block on the
      y_z-to-b_k block chain of Z) == both-parity simple z-to-b_k paths in
      B-(P-z), the latter by exhaustive DFS (step-capped).
  (2) GATES/BOUNDARY: every component X of the certificate zone S has exactly one
      crown-gate (0 in the degenerate S=Z case), and |N_B(X)| >= 3.
  (3) CENSUS: stratum of each configuration -- T3'-witnessed (feet of both
      parities), F2 monochrome (>=2 one parity), F1 (one foot), F0 (none);
      pigeonhole: in F0, S has one component and it takes both endpoint
      off-path edges; in F1, S has at most two components.
  (4) REPLICATION: each certificate z yields an even {i,k}- and an even
      {j,k}-connecting path by prefix/suffix + both-parity continuation.
  (5) REROUTE: for same-parity feet pairs z,z' of one S-component, the path
      P'' = P[b_i..z] + zy + W + y'z' + P[z'..b_j] with W an even y-to-y' path
      in Z is a simple connecting path (the Shadow Pin's carrier).
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0])

def tau_ge2(a):
    def bip(b):
        c={}
        for s in b:
            if s in c: continue
            c[s]=0; q=deque([s])
            while q:
                u=q.popleft()
                for w in b[u]:
                    if w not in c: c[w]=1-c[u]; q.append(w)
                    elif c[w]==c[u]: return False
        return True
    if bip(a): return False
    for e in {frozenset((u,w)) for u in a for w in a[u]}:
        u,v=tuple(e)
        b={x:set(a[x]) for x in a}; b[u].discard(v); b[v].discard(u)
        if bip(b): return False
    return True

def rand_conn(a,terms,i,j,rng):
    k=3-i-j; s,t=terms[i],terms[j]
    p=[s]; on={s,terms[k]}
    while p[-1]!=t:
        nb=[w for w in a[p[-1]] if w not in on]
        if not nb: return None
        w=rng.choice(nb); p.append(w); on.add(w)
        if len(p)>len(a): return None
    return p

def comp_of(a,verts,s):
    if s not in verts: return None
    seen={s}; q=deque([s])
    while q:
        u=q.popleft()
        for w in a[u]:
            if w in verts and w not in seen: seen.add(w); q.append(w)
    return seen

def induced(a,verts):
    return {v:{w for w in a[v] if w in verts} for v in verts}

def nonbip(g,verts):
    c={}
    for s in verts:
        if s in c: continue
        c[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in g[u]:
                if w not in c: c[w]=1-c[u]; q.append(w)
                elif c[w]==c[u]: return True
    return False

def blocks_of(g):
    """biconnected components (vertex sets) of graph dict g; iterative DFS."""
    idx={}; low={}; stack=[]; blocks=[]; t=[0]
    for root in g:
        if root in idx: continue
        st=[(root,None,iter(g[root]))]; idx[root]=low[root]=t[0]; t[0]+=1
        while st:
            v,parent,it=st[-1]
            adv=False
            for w in it:
                if w==parent: continue
                if w not in idx:
                    stack.append((v,w)); idx[w]=low[w]=t[0]; t[0]+=1
                    st.append((w,v,iter(g[w]))); adv=True; break
                elif idx[w]<idx[v]:
                    stack.append((v,w)); low[v]=min(low[v],idx[w])
            if adv: continue
            st.pop()
            if st:
                pv=st[-1][0]; low[pv]=min(low[pv],low[v])
                if low[v]>=idx[pv]:
                    blk=set();
                    while stack:
                        e=stack.pop(); blk.add(e[0]); blk.add(e[1])
                        if e==(pv,v): break
                    if blk: blocks.append(frozenset(blk))
    return blocks

def chain_blocks(g,blocks,y,target):
    """blocks on the block-cut tree path from y to target (inclusive)."""
    from collections import defaultdict
    inb=defaultdict(list)
    for bi,B in enumerate(blocks):
        for v in B: inb[v].append(bi)
    cuts={v for v in inb if len(inb[v])>1}
    T=defaultdict(set)
    for bi,B in enumerate(blocks):
        for v in B&cuts:
            T[('b',bi)].add(('c',v)); T[('c',v)].add(('b',bi))
    def node(v):
        if v in cuts: return ('c',v)
        return ('b',inb[v][0])
    s,tt=node(y),node(target)
    if s==tt:
        return [blocks[s[1]]] if s[0]=='b' else []
    par={s:None}; q=deque([s])
    while q:
        u=q.popleft()
        if u==tt: break
        for w in T[u]:
            if w not in par: par[w]=u; q.append(w)
    if tt not in par: return None
    path=[]; u=tt
    while u is not None: path.append(u); u=par[u]
    return [blocks[x[1]] for x in path if x[0]=='b']

def simple_path_parities(g,s,t,cap=250000):
    """parities of simple s-t paths in graph dict g, exhaustive DFS with step cap.
    returns (set_of_parities, complete_bool)"""
    par=set(); steps=[0]
    onpath={s}
    def dfs(u,d):
        if steps[0]>cap: return
        if u==t:
            par.add(d%2); return
        for w in g[u]:
            if w in onpath: continue
            steps[0]+=1
            onpath.add(w); dfs(w,d+1); onpath.discard(w)
            if len(par)==2: return
    dfs(s,0)
    return par, steps[0]<=cap

def rand_simple_path(g,s,t,rng,want_par=None,tries=1500):
    for _ in range(tries):
        p=[s]; on={s}
        while p[-1]!=t:
            nb=[w for w in g[p[-1]] if w not in on]
            if not nb: break
            w=rng.choice(nb); p.append(w); on.add(w)
        if p[-1]==t and (want_par is None or (len(p)-1)%2==want_par):
            return p
    return None

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 408)
made=0
tot={'cfg':0,'crit_ok':0,'crit_bad':0,'crit_skip':0,'gate_ok':0,'gate_bad':0,
     'bd_ok':0,'bd_bad':0,'T3':0,'F2':0,'F1':0,'F0':0,'F0_pig_ok':0,'F0_pig_bad':0,
     'F1_pig_ok':0,'F1_pig_bad':0,'rep_ok':0,'rep_fail':0,'rr_ok':0,'rr_fail':0,
     'cert2':0,'pairs':0}
while made<8:
    a,terms=make_real_bridge(30,rng)
    if a is None or not tau_ge2(a): continue
    made+=1
    line=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        k=3-i-j; bk=terms[k]; got=set()
        for _ in range(4000):
            if len(got)>=20: break
            P=rand_conn(a,terms,i,j,rng)
            if P is None or (len(P)-1)%2: continue
            key=tuple(P)
            if key in got or tuple(reversed(P)) in got: continue
            got.add(key)
        parities_seen=set()
        for key in got:
            P=list(key); Pset=set(P)
            Z=comp_of(a,set(a)-Pset,bk)
            if Z is None: continue
            gZ=induced(a,Z)
            if not nonbip(gZ,Z): continue
            tot['cfg']+=1
            blocks=blocks_of(gZ)
            oddb={B for B in blocks if nonbip(induced(a,B),B)}
            # certificate zone S and per-vertex chain criterion
            S=set()
            chain_ok={}
            for y in Z:
                ch=chain_blocks(gZ,blocks,y,bk)
                c=(ch is not None) and any(B in oddb for B in ch)
                chain_ok[y]=c
                if c: S.add(y)
            # (1) criterion check on interior vertices
            interior=P[1:-1]
            feet={}
            for idx,z in enumerate(interior):
                third=[w for w in a[z] if w!=P[idx] and w!=P[idx+2]]
                if len(third)!=1: continue
                y=third[0]
                if y not in Z: continue
                feet[z]=(idx+1,y)   # position of z along P, its Z-neighbour
                gz={v:set(gZ[v]) for v in Z}
                gz[z]={y}; gz[y]=gZ[y]|{z}
                par,complete=simple_path_parities(gz,z,bk)
                if not complete:
                    tot['crit_skip']+=1
                else:
                    both=(len(par)==2)
                    if both==chain_ok[y]: tot['crit_ok']+=1
                    else: tot['crit_bad']+=1
            # (2) gates and boundary of S-components
            Scomps=[]
            left=set(S)
            while left:
                s0=next(iter(left)); X=comp_of(gZ,S,s0)&S
                X=comp_of(induced(gZ,S),S,s0)
                Scomps.append(X); left-=X
            for X in Scomps:
                gates={w for v in X for w in gZ[v] if w not in X}
                nb=set()
                for v in X:
                    nb|={w for w in a[v] if w not in X}
                if bk in X:
                    ok_g=(len(gates)==0)
                else:
                    ok_g=(len(gates)==1)
                tot['gate_ok' if ok_g else 'gate_bad']+=1
                tot['bd_ok' if len(nb)>=3 else 'bd_bad']+=1
            # (3) census
            cert=[z for z in feet if chain_ok[feet[z][1]]]
            pars={feet[z][0]%2 for z in cert}
            parities_seen|=pars
            if len(pars)==2: tot['T3']+=1
            elif len(cert)>=2: tot['F2']+=1
            elif len(cert)==1:
                tot['F1']+=1
                tot['F1_pig_ok' if len(Scomps)<=2 else 'F1_pig_bad']+=1
            else:
                tot['F0']+=1
                # pigeonhole: one S-component, taking both endpoint off-edges
                bi,bj=terms[i],terms[j]
                offi=[w for w in a[bi] if w!=P[1]]; offj=[w for w in a[bj] if w!=P[-2]]
                ok=len(Scomps)==1 and offi and offj and offi[0] in Scomps[0] and offj[0] in Scomps[0]
                tot['F0_pig_ok' if ok else 'F0_pig_bad']+=1
            # (4) replication from up to 2 certificates
            for z in cert[:2]:
                pos,y=feet[z]
                gz={v:set(gZ[v]) for v in Z}
                gz[z]={y}; gz[y]=gZ[y]|{z}
                need=pos%2   # continuation parity making prefix+cont even
                W=rand_simple_path(gz,z,bk,rng,want_par=need)
                need2=(len(P)-1-pos)%2
                W2=rand_simple_path(gz,z,bk,rng,want_par=need2)
                if W and W2: tot['rep_ok']+=1
                else: tot['rep_fail']+=1
            # (5) reroute carrier for same-parity feet pairs in one S-component
            done=False
            for z1,z2 in itertools.combinations(cert,2):
                if done: break
                p1,y1=feet[z1]; p2,y2=feet[z2]
                if p1%2!=p2%2: continue
                X=next((X for X in Scomps if y1 in X),None)
                if X is None or y2 not in X: continue
                W=rand_simple_path(induced(gZ,X),y1,y2,rng,want_par=0)
                if W is None: continue
                if p1>p2: (p1,y1,z1),(p2,y2,z2)=(p2,y2,z2),(p1,y1,z1); W=W[::-1]
                Pp=P[:p1+1]+W+P[p2:]
                ok=(len(set(Pp))==len(Pp)) and (len(Pp)-1)%2==0 and Pp[0]==P[0] and Pp[-1]==P[-1]
                edges_ok=all(Pp[t+1] in a[Pp[t]] for t in range(len(Pp)-1))
                tot['rr_ok' if (ok and edges_ok) else 'rr_fail']+=1
                done=True
        tot['pairs']+=1
        if len(parities_seen)==2: tot['cert2']+=1
        line.append('%s'%('BOTH' if len(parities_seen)==2 else ('one' if parities_seen else 'none')))
    print('  bridge %d  cert-parities per pair: %s'%(made,'  '.join(line)))

print()
print('configs (Z nonbip): %d'%tot['cfg'])
print('(1) criterion:  ok %d   BAD %d   skipped(cap) %d'%(tot['crit_ok'],tot['crit_bad'],tot['crit_skip']))
print('(2) gates:      ok %d   BAD %d   |  boundary>=3: ok %d  BAD %d'%(tot['gate_ok'],tot['gate_bad'],tot['bd_ok'],tot['bd_bad']))
print('(3) census:     T3-witnessed %d   F2 %d   F1 %d   F0 %d'%(tot['T3'],tot['F2'],tot['F1'],tot['F0']))
print('    pigeonhole: F0 ok %d BAD %d   F1 ok %d BAD %d'%(tot['F0_pig_ok'],tot['F0_pig_bad'],tot['F1_pig_ok'],tot['F1_pig_bad']))
print('(4) replication: ok %d   fail %d'%(tot['rep_ok'],tot['rep_fail']))
print('(5) reroute:     ok %d   fail %d'%(tot['rr_ok'],tot['rr_fail']))
print('(CERT2) both parities seen: %d of %d pairs'%(tot['cert2'],tot['pairs']))
