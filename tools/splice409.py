"""splice409.py -- the certificate game's laws, audited.

(1) SPLICE INVARIANCE: single even-preserving splices of even connecting paths
    never flip any surviving position parity; parity-flipping splices always
    produce odd paths.
(2) FLIPPER: split-vertex witnesses (z certified on two even paths at opposite
    prefix parities) always have a non-bipartite prefix-union (the flipper odd
    cycle); count the shared-certifier pattern (certifying block avoids BOTH
    paths) vs the crossing pattern.
(3) PAIR SERVICE: an odd cycle through b_i avoiding b_k (b_i in a non-bipartite
    block of B - b_k) implies connecting paths of BOTH parities for {i,j}.
(4) SATURATION: if b_k's block in Z(P) is non-bipartite, EVERY interior foot
    into Z(P) is a certificate.
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
                    blk=set()
                    while stack:
                        e=stack.pop(); blk.add(e[0]); blk.add(e[1])
                        if e==(pv,v): break
                    if blk: blocks.append(frozenset(blk))
    return blocks

def chain_blocks(g,blocks,y,target):
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

def config(a,terms,P,bk):
    """returns (Z, certificates {z:(pos,y)}, oddblocks-on-chain map, blocks, sat_flag)"""
    Pset=set(P); Z=comp_of(a,set(a)-Pset,bk)
    if Z is None: return None
    gZ=induced(a,Z)
    if not nonbip(gZ,Z): return None
    blocks=blocks_of(gZ)
    oddb=[B for B in blocks if nonbip(induced(a,B),B)]
    oddset=set(oddb)
    cert={}; feet={}
    for idx in range(1,len(P)-1):
        z=P[idx]
        third=[w for w in a[z] if w!=P[idx-1] and w!=P[idx+1]]
        if len(third)!=1: continue
        y=third[0]
        if y not in Z: continue
        feet[z]=(idx,y)
        ch=chain_blocks(gZ,blocks,y,bk)
        if ch is not None and any(B in oddset for B in ch):
            qb=[B for B in ch if B in oddset]
            cert[z]=(idx,y,qb)
    bkodd=any(bk in B for B in oddb)
    return Z,feet,cert,bkodd

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 409)
made=0
T={'spl_even':0,'spl_flip':0,'spl_oddpath':0,'sv':0,'sv_flip_ok':0,'sv_flip_bad':0,
   'sv_shared':0,'sv_cross':0,'serv_hyp':0,'serv_ok':0,'serv_bad':0,
   'sat_cfg':0,'sat_ok':0,'sat_bad':0}
while made<8:
    a,terms=make_real_bridge(30,rng)
    if a is None or not tau_ge2(a): continue
    made+=1
    for i,j in [(0,1),(0,2),(1,2)]:
        k=3-i-j; bk=terms[k]
        # (3) pair service: hypothesis and conclusion
        rest=set(a)-{bk}
        gR=induced(a,rest)
        blR=blocks_of(gR)
        hyp=any(terms[i] in B and nonbip(induced(a,B),B) for B in blR)
        if hyp:
            T['serv_hyp']+=1
            pars=set()
            for _ in range(4000):
                P=rand_conn(a,terms,i,j,rng)
                if P: pars.add((len(P)-1)%2)
                if len(pars)==2: break
            T['serv_ok' if len(pars)==2 else 'serv_bad']+=1
        # sample even paths and configs
        paths=[]; seen=set()
        for _ in range(4000):
            if len(paths)>=15: break
            P=rand_conn(a,terms,i,j,rng)
            if P is None or (len(P)-1)%2: continue
            key=tuple(P)
            if key in seen or tuple(reversed(P)) in seen: continue
            seen.add(key); paths.append(P)
        cfgs=[]
        for P in paths:
            c=config(a,terms,P,bk)
            if c: cfgs.append((P,)+c)
        # (4) saturation
        for P,Z,feet,cert,bkodd in cfgs:
            if bkodd:
                T['sat_cfg']+=1
                ok=all(z in cert for z in feet)
                T['sat_ok' if ok else 'sat_bad']+=1
        # (1) splice invariance: reroute through a non-Z component between two feet
        done=0
        for P,Z,feet,cert,bkodd in cfgs:
            if done>=3: break
            Pset=set(P)
            comps=[]
            left=set(a)-Pset-Z
            while left:
                s0=next(iter(left)); D=comp_of(a,left|{s0},s0); comps.append(D); left-=D
            for D in comps:
                dfeet=[(idx,P[idx],w) for idx in range(1,len(P)-1) for w in a[P[idx]]
                       if w in D and P[idx+1]!=w and P[idx-1]!=w]
                if len(dfeet)<2: continue
                (m,zm,x),(m2,z2,x2)=dfeet[0],dfeet[-1]
                if m>=m2: continue
                gD=induced(a,D)
                W=None
                for _ in range(400):
                    p=[x]; on={x}
                    while p[-1]!=x2:
                        nb=[w for w in gD[p[-1]] if w not in on]
                        if not nb: break
                        w=rng.choice(nb); p.append(w); on.add(w)
                    if p[-1]==x2: W=p; break
                if W is None: continue
                ell=len(W)+1  # edges: zm-x, W path edges, x2-z2 => len(W)-1+2
                Pp=P[:m+1]+W+P[m2:]
                if len(set(Pp))!=len(Pp): continue
                newlen=len(Pp)-1; oldlen=len(P)-1
                if newlen%2==0:
                    T['spl_even']+=1
                    flip=False
                    for q in range(m2,len(P)):
                        oldpos=q; newpos=q-(m2-m)+(len(W)+1)
                        if (oldpos-newpos)%2: flip=True
                    if flip: T['spl_flip']+=1
                else:
                    T['spl_oddpath']+=1
                done+=1
                break
        # (2) split-vertex witnesses
        byv={}
        for ci,(P,Z,feet,cert,bkodd) in enumerate(cfgs):
            for z,(idx,y,qb) in cert.items():
                byv.setdefault(z,{}).setdefault(idx%2,(ci,idx))
        for z,d in byv.items():
            if len(d)<2: continue
            T['sv']+=1
            (c1,d1),(c2,d2)=d[0],d[1]
            P1=cfgs[c1][0]; P2=cfgs[c2][0]
            pre1=P1[:d1+1]; pre2=P2[:d2+1]
            U={}
            for pr in (pre1,pre2):
                for t in range(len(pr)-1):
                    U.setdefault(pr[t],set()).add(pr[t+1]); U.setdefault(pr[t+1],set()).add(pr[t])
            T['sv_flip_ok' if nonbip(U,set(U)) else 'sv_flip_bad']+=1
            # shared certifier: a certifying block of z on P1 avoiding P2 and inside Z(P2)
            qb1=cfgs[c1][3][z][2]; Z2=cfgs[c2][1]; P2set=set(P2)
            shared=any((B&P2set)==set() and B<=Z2 for B in qb1)
            T['sv_shared' if shared else 'sv_cross']+=1
    print('  bridge %d done'%made)

print()
print('(1) splice: even splices %d, parity flips %d (expect 0); parity-flip attempts -> odd path %d'
      %(T['spl_even'],T['spl_flip'],T['spl_oddpath']))
print('(2) split-vertex witnesses %d: flipper found %d, missing %d (expect 0 missing);'
      %(T['sv'],T['sv_flip_ok'],T['sv_flip_bad']))
print('    shared-certifier pattern %d, crossing pattern %d'%(T['sv_shared'],T['sv_cross']))
print('(3) pair service: hypothesis held %d, both parities found %d, failed %d (expect 0)'
      %(T['serv_hyp'],T['serv_ok'],T['serv_bad']))
print('(4) saturation: b_k-odd-block configs %d, all feet certified %d, violations %d (expect 0)'
      %(T['sat_cfg'],T['sat_ok'],T['sat_bad']))
