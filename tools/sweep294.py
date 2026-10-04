"""sweep294.py -- aggregate anatomy of (L2**)-failing pairs over many bridges."""
import sys, random, itertools, collections
src = open('tools/fail294.py').read()
src = src.split('def main()')[0].split('"""', 2)[2]
exec(src)

seeds = [int(x) for x in sys.argv[1].split(',')]
per   = int(sys.argv[2])
agg = {'graphs':0,'pairs':0,'ok':0,'fail':0,'recs':[]}
for s in seeds:
    rng = random.Random(s)
    for n in (13,15,17,19,21):
        got = 0
        while got < per:
            a, terms = random_bridge(n, rng)
            if a is None: break
            got += 1; agg['graphs'] += 1
            analyse(a, terms, agg)
R = agg['recs']
print('bridges %d   terminal pairs %d   (L2**) holds %d   FAILS %d'
      % (agg['graphs'], agg['pairs'], agg['ok'], agg['fail']))
if not R: sys.exit()
def tally(key):
    c = collections.Counter(str(r[key]) for r in R)
    return '  '.join('%s:%d' % kv for kv in sorted(c.items()))
print()
print('FAILING PAIRS, %d in all' % len(R))
for k in ['core_size','core_distance','core_separates','core_meets_all_odd',
          'nu_odd','split','SC','OCT_all_even_paths']:
    print('  %-22s %s' % (k, tally(k)))
tot_link  = sum(r['both_terminals_on_link'] for r in R)
tot_odd   = sum(r['both_arcs_odd'] for r in R)
tot_1link = sum(r['single_link'] for r in R)
tot_disc  = sum(r['even_paths_disconnecting'] for r in R)
tot_even  = sum(r['even'] for r in R)
print('  %-22s %d' % ('even conn. paths', tot_even))
print('  %-22s min %d  median %d  max %d'
      % ('per failing pair', min(r['even'] for r in R),
         sorted(r['even'] for r in R)[len(R)//2], max(r['even'] for r in R)))
print('  %-22s %d of %d even paths (250.2 says nothing about these)'
      % ('int(alpha) disconnects', tot_disc, tot_even))
print('  %-22s %d' % ('pairs P,R with 1 link', tot_1link))
print('  %-22s %d' % ('  .. b_1,b_2 on it', tot_link))
print('  %-22s %d   (both arcs EVEN would be (L2**))' % ('  .. both arcs ODD', tot_odd))
print('  %-22s %d' % ('second even path found', sum(r['second_even_path_exists'] for r in R)))
