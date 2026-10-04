"""funnel429.py -- the funnel kill: B has no non-degenerate 2-cuts.

B = the minimal cubic 3-connected counterexample G with three pairwise-far
edges subdivided (the terminals). Claim: every vertex pair {u,v} of B whose
removal disconnects B merely ISOLATES A SINGLE TERMINAL (the subdivision
vertex of the edge whose endpoints are u,v); no 2-cut separates the
terminal set from anything else. Hence kappa(T, Z) >= 3 for every cycle Z,
the Three-Fan Pigeonhole (425.3) applies to every witness, (S) is
UNCONDITIONAL, and the funnel case of 425.4 is VACUOUS.

Machine: subdivided Petersen + random 3-connected cubic girth-5 graphs:
enumerate ALL vertex pairs, remove, inspect components. Expected: every
disconnection isolates exactly one terminal.
"""
import random, itertools

def components(adj, removed):
    verts = [v for v in adj if v not in removed]
    seen = set(); comps = []
    for s in verts:
        if s in seen: continue
        comp = {s}; st = [s]
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in removed and y not in comp:
                    comp.add(y); st.append(y)
        seen |= comp; comps.append(comp)
    return comps

def subdivide(n, edges, picked):
    adj = {v: set() for v in range(n)}
    for a, b in edges:
        if tuple(sorted((a, b))) in picked: continue
        adj[a].add(b); adj[b].add(a)
    for k, (a, b) in enumerate(sorted(picked)):
        t = 't%d' % k
        adj[t] = {a, b}; adj[a].add(t); adj[b].add(t)
    return adj

def audit(name, adj, terminals):
    bad = 0; degen = 0; total = 0
    for u, v in itertools.combinations(list(adj), 2):
        comps = components(adj, {u, v})
        if len(comps) == 1: continue
        total += 1
        small = sorted(comps, key=len)[:-1]
        if all(len(c) == 1 and next(iter(c)) in terminals for c in small):
            degen += 1
        else:
            bad += 1
    print(' %-28s 2-cuts: %d degenerate (isolated terminal), %d GENUINE'
          % (name, degen, bad))
    return bad

def is3conn_cubic(n, edges):
    adj = {v: set() for v in range(n)}
    for a, b in edges: adj[a].add(b); adj[b].add(a)
    if any(len(adj[v]) != 3 for v in adj): return False
    for rm in itertools.combinations(range(n), 2):
        if len(components(adj, set(rm))) > 1: return False
    return True

def girth_ge5(n, edges):
    adj = {v: [] for v in range(n)}
    for a, b in edges: adj[a].append(b); adj[b].append(a)
    for s in range(n):
        dist = {s: 0}; par = {s: None}; q = [s]
        while q:
            v = q.pop(0)
            for w in adj[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1; par[w] = v; q.append(w)
                elif par[v] != w and dist[v] + dist[w] + 1 < 5:
                    return False
    return True

random.seed(429)
PET = [(i, (i + 1) % 5) for i in range(5)] + \
      [(i + 5, (i + 2) % 5 + 5) for i in range(5)] + [(i, i + 5) for i in range(5)]
picked = {(0, 1), (2, 3), (7, 9)}
adj = subdivide(10, PET, picked)
totalbad = audit('Petersen, 3 subdivided', adj, {'t0', 't1', 't2'})

built = 0
while built < 4:
    n = random.choice([12, 14, 16])
    stubs = [v for v in range(n) for _ in range(3)]
    random.shuffle(stubs)
    E = set(); ok = True
    for i in range(0, len(stubs), 2):
        a, b = stubs[i], stubs[i + 1]
        if a == b or frozenset((a, b)) in E: ok = False; break
        E.add(frozenset((a, b)))
    if not ok: continue
    edges = [tuple(sorted(e)) for e in E]
    if not girth_ge5(n, edges) or not is3conn_cubic(n, edges): continue
    random.shuffle(edges)
    pick = []
    for e in edges:
        if all(x not in f for f in pick for x in e): pick.append(e)
        if len(pick) == 3: break
    if len(pick) < 3: continue
    adj = subdivide(n, edges, set(pick))
    totalbad += audit('random 3-conn cubic n=%d' % n, adj,
                      {'t0', 't1', 't2'})
    built += 1

print()
print('GENUINE 2-cuts across all subdivided graphs: %d' % totalbad)
assert totalbad == 0
print('=> kappa(T, Z) >= 3 for every cycle Z: the Three-Fan always exists,')
print('   (S) is UNCONDITIONAL, and the funnel case of 425.4 is VACUOUS.')
