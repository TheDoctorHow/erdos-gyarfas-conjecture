"""anatomy.py -- dissect the smallest (L2**)-failing pairs: where is the critical edge,
and what does B - e look like?"""
import sys, random, itertools
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split('"""', 2)[2])

def all_paths(a, s, t):
    out = []; path = [s]; onp = {s}
    def rec(u):
        for w in a[u]:
            if w == t: out.append(path + [t])
            elif w not in onp:
                onp.add(w); path.append(w); rec(w); path.pop(); onp.discard(w)
    rec(s); return out

def blocks(a):
    """blocks (2-connected components) of a, as edge sets, plus cut vertices"""
    idx = {}; low = {}; par = {}; st = []; out = []; cuts = set(); C = [0]
    def dfs(u):
        idx[u] = low[u] = C[0]; C[0] += 1; ch = 0
        for w in a[u]:
            if w not in idx:
                ch += 1; par[w] = u; st.append(frozenset((u, w))); dfs(w)
                low[u] = min(low[u], low[w])
                if low[w] >= idx[u]:
                    if par.get(u) is not None or ch > 1: cuts.add(u)
                    blk = []
                    while True:
                        e = st.pop(); blk.append(e)
                        if e == frozenset((u, w)): break
                    out.append(blk)
            elif w != par.get(u) and idx[w] < idx[u]:
                low[u] = min(low[u], idx[w]); st.append(frozenset((u, w)))
    for s in a:
        if s not in idx: par[s] = None; dfs(s)
    return out, cuts

def bip(es):
    a = adj_of(es); col = {}
    for s in a:
        if s in col: continue
        col[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in col: col[w] = 1 - col[u]; q.append(w)
                elif col[w] == col[u]: return False
    return True

def report(a, terms, b1, b2, ev):
    core = set(eset(ev[0]))
    for p in ev[1:]: core &= eset(p)
    E = {frozenset((u, w)) for u in a for w in a[u]}
    print('n=%d  m=%d  terminals=%s  pair=(%s,%s)  #even=%d'
          % (len(a), len(E), terms, b1, b2, len(ev)))
    print('  edges:', ' '.join('%s%s' % tuple(sorted(map(str, e))) for e in
                               sorted(E, key=lambda s: sorted(map(str, s)))))
    for e in core:
        u, v = tuple(e)
        print('  CRITICAL EDGE  %s-%s' % (u, v))
        H = E - {e}
        aH = adj_of(H)
        print('    B-e bipartite overall :', bip(H))
        bl, cuts = blocks(aH)
        print('    B-e has %d blocks, cut vertices %s'
              % (len(bl), sorted(map(str, cuts))))
        for k, blk in enumerate(sorted(bl, key=len)):
            vs = {x for f in blk for x in f}
            marks = [str(z) for z in (b1, b2, u, v) if z in vs]
            t3 = [str(z) for z in terms if z in vs]
            print('      block %d: %2d edges  bipartite=%-5s  contains %s'
                  % (k, len(blk), bip(blk), ' '.join(sorted(set(marks + t3))) or '-'))
        # all odd cycles: do they all use e?
        odd_avoiding = 0
        for x in aH:
            pass
        print('    every u-v path in B-e even? ', end='')
        d = {(u, 0)}; q = deque([(u, 0)])
        while q:
            z, p = q.popleft()
            for w in aH.get(z, ()):
                if (w, 1 - p) not in d: d.add((w, 1 - p)); q.append((w, 1 - p))
        print('even-reach=%s odd-reach=%s' % ((v, 0) in d, (v, 1) in d))
        print('    b_1-b_2 in B-e:  even-reach/odd-reach = ', end='')
        d = {(b1, 0)}; q = deque([(b1, 0)])
        while q:
            z, p = q.popleft()
            for w in aH.get(z, ()):
                if (w, 1 - p) not in d: d.add((w, 1 - p)); q.append((w, 1 - p))
        print((b2, 0) in d, (b2, 1) in d)

def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 2)
    want = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    shown = 0
    for m in (10, 12, 14):
        for _ in range(4000):
            if shown >= want: return
            a, terms = make_bridge(m, rng)
            if a is None: continue
            for b1, b2 in itertools.combinations(terms, 2):
                ps = all_paths(a, b1, b2)
                ev = [p for p in ps if (len(p) - 1) % 2 == 0]
                if len(ev) < 2: continue
                core = set(eset(ev[0]))
                for p in ev[1:]: core &= eset(p)
                if not core: continue
                print('=' * 72); report(a, terms, b1, b2, ev); shown += 1
                if shown >= want: return

main()
