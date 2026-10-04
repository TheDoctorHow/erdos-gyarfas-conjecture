import random, itertools, sys

def girth(adj):
    n = len(adj); best = 10**9
    for s in range(n):
        dist = [-1]*n; par = [-1]*n; dist[s]=0; q=[s]
        for u in q:
            for v in adj[u]:
                if dist[v] < 0:
                    dist[v]=dist[u]+1; par[v]=u; q.append(v)
                elif v != par[u]:
                    best = min(best, dist[u]+dist[v]+1)
    return best

def connected_after(adj, removed):
    n = len(adj)
    start = next((v for v in range(n) if v not in removed), None)
    if start is None: return True
    seen = {start}; st=[start]
    while st:
        u = st.pop()
        for v in adj[u]:
            if v not in removed and v not in seen:
                seen.add(v); st.append(v)
    return len(seen) == n - len(removed)

def three_connected(adj):
    n = len(adj)
    if n < 4: return False
    for pair in itertools.combinations(range(n), 2):
        if not connected_after(adj, set(pair)): return False
    return True

def count_cycles(adj):
    """Return dict length -> number of cycles."""
    n = len(adj); res = {}
    for s in range(n):
        # paths starting at s, all other vertices > s
        stack = [(s, [s], {s})]
        while stack:
            u, path, seen = stack.pop()
            for v in adj[u]:
                if v == s and len(path) >= 3:
                    res[len(path)] = res.get(len(path), 0) + 1
                elif v > s and v not in seen:
                    stack.append((v, path+[v], seen | {v}))
    return {L: c//2 for L, c in res.items()}

def report(name, adj):
    g = girth(adj); tc = three_connected(adj)
    cyc = count_cycles(adj)
    N = [0]*8
    for L, c in cyc.items(): N[L % 8] += c
    total = sum(cyc.values())
    print("%-22s n=%-3d girth=%-2d 3conn=%-5s N0=%-6d N0mod2=%d total=%-7d N=%s"
          % (name, len(adj), g, tc, N[0], N[0] % 2, total, N))
    return N

def from_edges(n, edges):
    adj = [[] for _ in range(n)]
    for a,b in edges:
        adj[a].append(b); adj[b].append(a)
    return adj

def gp(n, k):
    """Generalized Petersen GP(n,k)."""
    E = []
    for i in range(n):
        E.append((i, (i+1) % n))
        E.append((i, n+i))
        E.append((n+i, n+((i+k) % n)))
    E = sorted({tuple(sorted(e)) for e in E})
    return from_edges(2*n, E)

def lcf(n, shifts):
    """LCF notation: Hamiltonian cycle 0..n-1 plus chords i -> i+shifts[i mod len]."""
    E = [(i, (i+1) % n) for i in range(n)]
    for i in range(n):
        s = shifts[i % len(shifts)]
        E.append((i, (i+s) % n))
    E = sorted({tuple(sorted(e)) for e in E})
    return from_edges(n, E)

print("=== named cubic graphs ===")
report("K4",           from_edges(4, [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]))
report("K33",          from_edges(6, [(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)]))
report("Q3(cube)",     lcf(8, [3,-3]))
report("Petersen",     gp(5,2))
report("Heawood",      lcf(14, [5,-5]))
report("MobiusKantor", lcf(16, [5,-5]))
report("Pappus",       lcf(18, [5,7,-7,7,-7,-5]))
report("Desargues",    lcf(20, [5,-5,9,-9]))
report("Dodecahedron", lcf(20, [10,7,4,-4,-7,10,-4,7,-7,4]))
report("Nauru",        lcf(24, [5,-9,7,-7,9,-5]))
report("McGee",        lcf(24, [12,7,-7]))

print()
print("=== generalized Petersen GP(n,2), girth>=5 ===")
for m in range(7, 13):
    report("GP(%d,2)" % m, gp(m,2))

print()
print("=== random cubic hamiltonian, girth>=5 ===")
random.seed(20260805)
found = 0
for n in range(10, 23, 2):
    got = 0
    tries = 0
    while got < 4 and tries < 4000:
        tries += 1
        perm = list(range(n)); random.shuffle(perm)
        pairs = [(perm[2*i], perm[2*i+1]) for i in range(n//2)]
        E = [(i,(i+1) % n) for i in range(n)]
        ok = True
        for a,b in pairs:
            if (b-a) % n in (1, n-1): ok = False; break
        if not ok: continue
        E = E + [tuple(sorted(p)) for p in pairs]
        if len(set(E)) != len(E): continue
        adj = from_edges(n, sorted(set(E)))
        if any(len(x) != 3 for x in adj): continue
        if girth(adj) < 5: continue
        if not three_connected(adj): continue
        got += 1; found += 1
        report("rand n=%d #%d" % (n, got), adj)
print("total random graphs tested:", found)
