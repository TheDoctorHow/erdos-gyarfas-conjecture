"""The extremal configuration at 2N = 18 is the SPARSE UNIFORM LADDER
S(2N, s, p) = { (pj + 1, pj + 1 + s) : j } -- chords of one span placed on
every p-th position -- together with extra chords inside the residue
class it leaves free.

This tests, for each (2N, s, p), whether the sparse ladder is girth-legal
and L-free, and then how many extra chords can be packed into the free
class, i.e. what margin the construction certifies.
"""
import sys, time, random
from pushL import short, has_len_through, bigpow

def sparse(n, s, p):
    ch = []
    used = set()
    j = 0
    while p * j + 1 < n + 1:
        a = (p * j + 1) % n
        b = (a + s) % n
        if a in used or b in used or a == b:
            break
        ch.append((min(a, b), max(a, b)))
        used.add(a); used.add(b)
        j += 1
        if p * j >= n:
            break
    return ch, used

def build(n, ch):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    for a, b in ch:
        adj[a].append(b); adj[b].append(a)
    return adj

def legal_and_free(n, ch, L):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    for a, b in ch:
        if short(adj, n, a, b):
            return None
        adj[a].append(b); adj[b].append(a)
        if has_len_through(adj, n, a, b, L):
            return None
    return adj

def fill(n, L, ch, rng, tries=200):
    """most extra chords addable on the free vertices"""
    best = list(ch)
    for _ in range(tries):
        cur = list(ch)
        adj = legal_and_free(n, cur, L)
        if adj is None:
            return None
        mate = [-1] * n
        for a, b in cur:
            mate[a] = b; mate[b] = a
        free = [v for v in range(n) if mate[v] == -1]
        rng.shuffle(free)
        for v in free:
            if mate[v] != -1:
                continue
            cand = [w for w in range(n) if w != v and mate[w] == -1
                    and not short(adj, n, v, w)]
            rng.shuffle(cand)
            for w in cand:
                adj[v].append(w); adj[w].append(v)
                if has_len_through(adj, n, v, w, L):
                    adj[v].pop(); adj[w].pop(); continue
                mate[v] = w; mate[w] = v
                cur.append((min(v, w), max(v, w)))
                break
        if len(cur) > len(best):
            best = cur
    return best

if __name__ == '__main__':
    rng = random.Random(5)
    for n in [int(x) for x in sys.argv[1:]]:
        L = bigpow(n)
        if L == n:
            print("  2N=%2d trivial" % n); continue
        N = n // 2
        bestrec = None
        for p in (2, 3, 4, 5):
            for s in range(4, n // 2 + 1):
                ch, _ = sparse(n, s, p)
                if len(ch) * 2 > n:
                    continue
                if legal_and_free(n, ch, L) is None:
                    continue
                full = fill(n, L, ch, rng)
                if full is None:
                    continue
                if bestrec is None or len(full) > bestrec[0]:
                    bestrec = (len(full), p, s, len(ch), full)
        if bestrec is None:
            print("  2N=%2d  L=%2d  no sparse uniform ladder works" % (n, L))
        else:
            k, p, s, base, full = bestrec
            print("  2N=%2d  N=%2d  L=%2d   S(2N,%d,p=%d) has %d chords, "
                  "filled to %d  =>  margin <= %d"
                  % (n, N, L, s, p, base, k, N - k))
            print("        ", sorted(full))
        sys.stdout.flush()
