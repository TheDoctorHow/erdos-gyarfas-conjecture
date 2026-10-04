"""fail294.py -- anatomy of genuine (L2**)-failing pairs.

Population: 214.2's model.  Take a 2-connected cubic graph of girth >= 5, subdivide
three edges; the subdivision points are the three terminals of degree 2, everything
else has degree 3.  Keep only non-bipartite, 2-connected, girth >= 5 results.

For every FAILING pair (no two internally disjoint even connecting paths) measure:
  * counts of even / odd connecting paths;
  * (UV)  is there a vertex internal to EVERY even connecting path?
  * (OCT) is the interior of every even connecting path an odd cycle transversal?
          -- 250.3/253.4(c) assert this; the deduction from 250.2 in fact needs
             b_1,b_2 to stay connected in B - int(alpha) with a non-bipartite block on
             their route, so this is a genuine test of an asserted fact;
  * whether b_1,b_2 remain connected in B - int(alpha);
  * nu_odd (max vertex-disjoint odd cycles) -- 254.1 predicts 1;
  * (SPLIT) do the even paths use both edges at each terminal?
  * link statistics: number of link cycles, whether b_1,b_2 share one, arc parities;
  * parity of |E(P) n E(R)| over even pairs, and whether (SC) holds for the pair.
"""
import sys, random, itertools
from collections import deque

exec(open('tools/links294.py').read().split('# ------------------------------------------------------------------ tests')[0]
     .split('"""', 2)[2])

def odd_cycles(a, cap=200000):
    """all odd cycles as vertex sets (small graphs only)"""
    vs = sorted(a, key=str); out = []
    for i, s in enumerate(vs):
        # cycles whose smallest vertex (by index) is s
        idx = {v: k for k, v in enumerate(vs)}
        path = [s]; onp = {s}
        def rec(u):
            if len(out) >= cap: return
            for w in a[u]:
                if idx[w] < i: continue
                if w == s and len(path) >= 3:
                    if len(path) % 2: out.append(frozenset(path))
                elif w not in onp:
                    onp.add(w); path.append(w); rec(w); path.pop(); onp.discard(w)
        rec(s)
    return list(set(out))

def nu_odd(oc):
    best = 0
    oc = sorted(oc, key=len)
    def rec(i, chosen):
        nonlocal best
        best = max(best, len(chosen))
        for j in range(i, len(oc)):
            if all(not (oc[j] & c) for c in chosen):
                rec(j + 1, chosen + [oc[j]])
    rec(0, [])
    return best

def blocks_route_nonbip(a, s, t, removed):
    """b_1,b_2 connected in a - removed?  and is some even s-t path present?"""
    sub = {v: {w for w in a[v] if w not in removed} for v in a if v not in removed}
    # BFS parity reachability
    seen = {(s, 0)}; q = deque([(s, 0)])
    while q:
        u, p = q.popleft()
        for w in sub.get(u, ()):
            if (w, 1 - p) not in seen:
                seen.add((w, 1 - p)); q.append((w, 1 - p))
    conn = any((t, p) in seen for p in (0, 1))
    even = (t, 0) in seen
    return conn, even

def analyse(a, terms, out):
    for b1, b2 in itertools.combinations(terms, 2):
        paths = all_paths(a, b1, b2)
        ev = [p for p in paths if (len(p) - 1) % 2 == 0]
        if not ev: continue
        fail = True
        for i in range(len(ev)):
            for j in range(i + 1, len(ev)):
                if not (eset(ev[i]) & eset(ev[j])):
                    fail = False; break
            if not fail: break
        out['pairs'] += 1
        if not fail:
            out['ok'] += 1
            continue
        out['fail'] += 1
        od = [p for p in paths if (len(p) - 1) % 2]
        rec = {'n': len(a), 'even': len(ev), 'odd': len(od),
               'len_even': sorted({len(p) - 1 for p in ev})}
        # universal vertex (common core of the even-path family)
        uni = set(ev[0][1:-1])
        for p in ev[1:]: uni &= set(p[1:-1])
        rec['core_size'] = len(uni)
        rec['universal_vertices'] = sorted(map(str, uni))
        # does the core also meet every ODD connecting path?
        rec['core_meets_all_odd'] = all(uni & set(p[1:-1]) for p in od) if od else None
        # is the core a b_1-b_2 separator?
        rec['core_separates'] = not blocks_route_nonbip(a, b1, b2, uni)[0] if uni else False
        # common core of ALL connecting paths
        allc = set(paths[0][1:-1])
        for p in paths[1:]: allc &= set(p[1:-1])
        rec['core_all_paths'] = sorted(map(str, allc))
        # are the core vertices adjacent / how far apart
        if len(uni) == 2:
            x, y = sorted(uni, key=str)
            rec['core_distance'] = dist(a, x, y)
        else:
            rec['core_distance'] = None
        # SPLIT
        rec['split'] = (len({p[1] for p in ev}) == 2 and len({p[-2] for p in ev}) == 2)
        # SC
        rec['SC'] = len({(len(p) - 1) % 4 for p in ev}) == 1
        # odd cycles / nu_odd / OCT
        oc = odd_cycles(a)
        rec['nu_odd'] = nu_odd(oc)
        oct_all = True; disc = 0; nonbip_route = 0
        for p in ev:
            interior = set(p[1:-1])
            if any(not (interior & c) for c in oc): oct_all = False
            conn, even2 = blocks_route_nonbip(a, b1, b2, interior)
            if not conn: disc += 1
            if even2: nonbip_route += 1
        rec['OCT_all_even_paths'] = oct_all
        rec['even_paths_disconnecting'] = disc
        rec['second_even_path_exists'] = nonbip_route      # must be 0 if failing
        # link statistics over even pairs
        one_link = both_on_link = both_arcs_odd = 0
        inter_odd = inter_even = 0
        for P, R in itertools.combinations(ev, 2):
            D = eset(P) ^ eset(R)
            cs = comps(D)
            c = len(eset(P) & eset(R))
            if c % 2: inter_odd += 1
            else: inter_even += 1
            if len(cs) == 1:
                one_link += 1
                blk, be = cs[0]
                if b1 in blk and b2 in blk:
                    both_on_link += 1
                    da = adj_of(be)
                    d = {b1: 0}; q = deque([b1])
                    while q:
                        u = q.popleft()
                        for w in da[u]:
                            if w not in d:
                                d[w] = d[u] + 1; q.append(w)
                    arc1 = d[b2]; arc2 = len(be) - arc1
                    if arc1 % 2 and arc2 % 2: both_arcs_odd += 1
        rec['even_pairs'] = one_link + (len(ev) * (len(ev) - 1) // 2 - one_link)
        rec['single_link'] = one_link
        rec['both_terminals_on_link'] = both_on_link
        rec['both_arcs_odd'] = both_arcs_odd
        rec['shared_edges_odd'] = inter_odd
        rec['shared_edges_even'] = inter_even
        out['recs'].append(rec)

def main():
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 20260805)
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    out = {'pairs': 0, 'ok': 0, 'fail': 0, 'graphs': 0, 'recs': []}
    for n in (13, 15, 17, 19):
        got = 0
        while got < per:
            a, terms = random_bridge(n, rng)
            if a is None: break
            got += 1; out['graphs'] += 1
            analyse(a, terms, out)
    print('graphs %d   pairs %d   (L2**) holds %d   FAILS %d'
          % (out['graphs'], out['pairs'], out['ok'], out['fail']))
    keys = ['n', 'even', 'odd', 'SC', 'split', 'nu_odd', 'core_size', 'core_distance', 'core_meets_all_odd', 'core_separates', 'core_all_paths',
            'OCT_all_even_paths', 'even_paths_disconnecting',
            'second_even_path_exists', 'single_link', 'both_terminals_on_link',
            'both_arcs_odd', 'shared_edges_odd', 'shared_edges_even']
    for r in out['recs']:
        print('-' * 60)
        for k in keys:
            print('  %-26s %s' % (k, r[k]))

main()
