"""budget309.py -- how far does the 308 engine reach as the degree-2 budget grows?

308 needed two numbers to meet:
  DESCENT THRESHOLD  T(s): if all handles are trees then every handle vertex has degree 2,
     so with s degree-2 vertices, two of them the pair a,b on the path P', one gets
     |V(K)| <= 3s - 6.  Hence |V(K)| >= 3s-5  ==>  some handle carries a cycle.
  BASE FLOOR  F(s): the least n admitting a degree sequence for a 2-connected bipartite
     subcubic graph of girth >= 6 with s vertices of degree 2 (cherry bound against girth).
The engine closes the stratum iff F(s) >= T(s), i.e. F(s) >= 3s-5.
"""
def floor_n(s, nmax=60):
    for n in range(3, nmax+1):
        for p in range(1, n):
            q = n - p
            for p2 in range(0, min(p, s)+1):
                q2 = s - p2
                if q2 < 0 or q2 > q: continue
                p3 = p - p2; q3 = q - q2
                if p3 < 0 or q3 < 0: continue
                if 3*p3 + 2*p2 != 3*q3 + 2*q2: continue
                m = 3*p3 + 2*p2
                if m < n: continue                       # 2-connected
                if 3*q3 + q2 > p*(p-1)//2: continue      # girth >= 6
                if 3*p3 + p2 > q*(q-1)//2: continue
                return n
    return None

print('   s   F(s)=floor   T(s)=3s-5   engine closes?   tau it covers')
for s in range(2, 16):
    F = floor_n(s); T = 3*s - 5
    tau = (s - 3) / 2.0
    print('  %2d    %-10s %-11d %-16s %s'
          % (s, F, T, 'YES' if (F is not None and F >= T) else 'no',
             ('tau <= %.1f' % tau) if s >= 3 else '-'))
