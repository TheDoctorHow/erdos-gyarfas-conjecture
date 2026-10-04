"""ecore294.py -- test the E-CORE conjecture.

  (E-CORE)  a terminal pair FAILS (L2**)  <=>  some EDGE lies on EVERY even
            connecting path between them.

The <= direction is trivial (two edge-disjoint even paths cannot both use e, and by
the Shared-Edge Lemma edge-disjoint = internally disjoint).  The => direction is a
HELLY property for the family of even connecting paths, and is what is tested here.
"""
import sys, random, itertools, collections
src = open('tools/fail294.py').read().split('def analyse')[0].split('"""', 2)[2]
exec(src)

seeds = [int(x) for x in sys.argv[1].split(',')]
per   = int(sys.argv[2])
ns    = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [13,15,17,19,21]

st = collections.Counter()
detail = []
for s in seeds:
    rng = random.Random(s)
    for n in ns:
        got = 0
        while got < per:
            a, terms = random_bridge(n, rng)
            if a is None: break
            got += 1; st['bridges'] += 1
            for b1, b2 in itertools.combinations(terms, 2):
                paths = all_paths(a, b1, b2)
                ev = [p for p in paths if (len(p)-1) % 2 == 0]
                if not ev: st['no_even_path'] += 1; continue
                st['pairs'] += 1
                fail = not any(not (eset(ev[i]) & eset(ev[j]))
                               for i in range(len(ev)) for j in range(i+1, len(ev)))
                # common edges of the whole even family
                core = set(eset(ev[0]))
                for p in ev[1:]: core &= eset(p)
                # common vertices
                vcore = set(ev[0][1:-1])
                for p in ev[1:]: vcore &= set(p[1:-1])
                if fail:
                    st['FAIL'] += 1
                    st['fail_has_edge_core' if core else 'fail_NO_edge_core'] += 1
                    if core:
                        st['edge_core_size_%d' % len(core)] += 1
                        # is the common edge incident to a terminal?
                        inc = sum(1 for e in core if b1 in e or b2 in e)
                        st['core_edges_at_terminal_%d' % inc] += 1
                    detail.append((n, len(ev), len(core), len(vcore)))
                else:
                    st['HOLDS'] += 1
                    if core: st['holds_has_edge_core'] += 1
print('bridges %d   pairs %d   (L2**) HOLDS %d   FAILS %d'
      % (st['bridges'], st['pairs'], st['HOLDS'], st['FAIL']))
for k in sorted(st):
    if k in ('bridges','pairs','HOLDS','FAIL'): continue
    print('  %-28s %d' % (k, st[k]))
if detail:
    print('\n  (n, #even paths, #common EDGES, #common VERTICES) for each failing pair:')
    for d in detail: print('   ', d)
