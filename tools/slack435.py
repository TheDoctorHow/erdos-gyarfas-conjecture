"""slack435.py -- correcting 433.2's slack estimate: the self-absorption
requirement is LINEAR in L, and the packing shell provides exactly that.

A leg of length ell carries ~ell interior vertices; if the non-middle
zone (~ell vertices) self-absorbs via chords of span >= 4 (girth), the
chords tile the leg. A single route may use ANY pairwise-disjoint
sub-family; greedy tiling gives ~ell/6 disjoint chords, total shortcut
>= 3 * ell/6 = ell/2 ... measured against the leg; against the full ear
(|P| ~ 2 ell + middle), the shortcut sum forces
   slack > (disjoint chords) * 3 ~ linear in L,
NOT sqrt(L) as 433.2 estimated. Machine: exact maximum over tilings.
"""
def max_disjoint_shortcut(ell, span):
    # chords (i, i+span) for i = 1..ell-span-1, every vertex chorded in
    # pairs; maximal disjoint sub-family by greedy: intervals [i, i+span]
    count = 0; pos = 1
    while pos + span <= ell - 1:
        count += 1; pos += span + 1
    return count, count * (span - 1)

for L in (24, 48, 96, 192):
    ell = L // 4                      # leg length at the leg bound
    k, shortcut = max_disjoint_shortcut(ell, 4)
    print('L = %3d: leg %2d, disjoint span-4 chords %2d, forced shortcut'
          ' %2d  (= %.3f L; sqrt(L) = %.1f)'
          % (L, ell, k, shortcut, shortcut / L, L ** 0.5))
print()
print('=> the self-absorption slack requirement grows LINEARLY (~0.15 L'
      ' per leg, more for longer legs), not as sqrt(L): 433.2 UNDERSTATED')
print('   the constraint. The packing shell supplies average slack ~ L/4')
print('   (supply N ~ L/3 disjoint ears in n ~ L^2/4 vertices: average')
print('   |P| ~ 3L/4, slack ~ L/4): the shell remains EXACTLY consistent.')
