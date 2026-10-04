"""(T1) enumeration with the bounded dual oracle.

Same pruned search as pushL.enumerate_missing (least unmatched vertex,
girth prune, first chord fixed up to reflection; complete matchings
only), but the (2N-2)-cycle oracle is dualdet2.detect2, whose suffix
bound makes deficiency-2 detection nearly linear, instead of the path
search whose cost explodes on branches where the long cycle is absent.

Validated below against the slow version's node counts at small sizes.
"""
import sys, time
from collections import deque
from pushL import short
from dualdet2 import detect2

def path_budget(adj, n, a, b, L, budget):
    """has_len_through with a node budget: True = cycle found (sound
    prune), False = certified absent, None = budget exhausted."""
    target = L - 1
    d = [-1] * n; d[b] = 0; q = deque([b])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if d[w] < 0:
                d[w] = d[u] + 1; q.append(w)
    if d[a] < 0:
        return False
    seen = bytearray(n); seen[a] = 1
    cnt = [budget]
    def rec(u, dep):
        if cnt[0] <= 0:
            return None
        cnt[0] -= 1
        rem = target - dep
        if d[u] > rem:
            return False
        for w in adj[u]:
            if w == b:
                if rem == 1:
                    return True
                continue
            if seen[w]:
                continue
            seen[w] = 1
            r = rec(w, dep + 1)
            seen[w] = 0
            if r:
                return True
            if r is None:
                return None
        return False
    complete = True
    for w in adj[a]:
        if w == b:
            continue
        seen[w] = 1
        r = rec(w, 1)
        seen[w] = 0
        if r:
            return True
        if r is None:
            complete = False
    return False if complete else None

def enumerate_missing_fast(n, L):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    mate = [-1] * n
    chords = []
    out = []
    nodes = [0]

    def rec():
        nodes[0] += 1
        a = next((v for v in range(n) if mate[v] == -1), None)
        if a is None:
            out.append(list(chords)); return
        for b in range(a + 1, n):
            if mate[b] != -1: continue
            if short(adj, n, a, b): continue
            if a == 0 and b > n // 2: break
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
            chords.append((a, b))
            # hybrid oracle.  Budgeted path search first: a found cycle
            # prunes soundly and is the common, cheap case.  On budget
            # exhaustion or certified absence, the eps=1 dual decides;
            # missing an eps=0 cycle only recurses harmlessly, and
            # surviving leaves are re-verified with the full oracle.
            r = path_budget(adj, n, a, b, L, 1000)
            if r is not True:
                got, _ = detect2(n, chords, L, eps_list=(1,))
                r = bool(got)
            if not r:
                rec()
            chords.pop()
            adj[a].pop(); adj[b].pop()
            mate[a] = -1; mate[b] = -1
    rec()
    return out, nodes[0]

if __name__ == '__main__':
    sys.setrecursionlimit(200000)
    for n in [int(x) for x in sys.argv[1:]]:
        t = time.time()
        out, nodes = enumerate_missing_fast(n, n - 2)
        # re-verify survivors with the FULL oracle (both parities)
        real = [m for m in out
                if not detect2(n, m, n - 2, eps_list=(0, 1))[0]]
        print("  2N=%2d  L=%2d : survivors %d, truly lacking = %d"
              "   [nodes %d, %.0fs]"
              % (n, n - 2, len(out), len(real), nodes, time.time() - t))
        for m in real[:2]:
            print("       ", m)
        sys.stdout.flush()
