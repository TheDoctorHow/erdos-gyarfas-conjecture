#!/bin/sh
# sweep.sh <n> <shards>   -- complete search of cubic C4-free graphs on n vertices
#   PART A: triangle-free  => girth >= 5 on n vertices
#   PART B: has a triangle => inflate(H,v) for H cubic on n-2 vertices
N=$1; S=$2
SH=sh/snarkhunter-2.0b
mkdir -p sw; rm -f sw/near_$N.txt sw/stat_$N.txt
i=0
while [ $i -lt $S ]; do
 ( cd $SH && ./snarkhunter.exe $N 5 o g s m $i $S 2>/dev/null ) \
   | EG_NEARFILE=sw/nearA_${N}_$i.txt ./egtool.exe -q >> sw/logA_${N}_$i.txt 2>> sw/statA_${N}_$i.txt &
 i=$((i+1))
done
i=0
while [ $i -lt $S ]; do
 ( cd $SH && ./snarkhunter.exe $((N-2)) 3 o g s m $i $S 2>/dev/null ) \
   | EG_NEARFILE=sw/nearB_${N}_$i.txt ./inflate1.exe >> sw/logB_${N}_$i.txt 2>> sw/statB_${N}_$i.txt &
 i=$((i+1))
done
wait
echo "=== n=$N COMPLETE ==="
echo "-- PART A (girth>=5) --"; grep -h "^graphs=" sw/statA_${N}_*.txt | awk -F'[= ]' '{g+=$2;nm+=$4;c+=$6} END{print "graphs="g" nearmiss="nm" counterexamples="c}'
echo "-- PART B (triangle) --"; grep -h "^INF1 " sw/statB_${N}_*.txt | awk '{for(i=2;i<=NF;i++){split($i,a,"=");v[a[1]]+=a[2]}} END{printf "bases=%d cand=%d noC4=%d noC4C8=%d ctr=%d\n",v["bases"],v["cand"],v["noC4"],v["noC4C8"],v["ctr"]}'
cat sw/nearA_${N}_*.txt sw/nearB_${N}_*.txt 2>/dev/null > sw/near_$N.txt
echo "-- near-miss records (may repeat isomorphic copies): $(wc -l < sw/near_$N.txt)"
grep -h COUNTEREXAMPLE sw/logA_${N}_*.txt sw/logB_${N}_*.txt 2>/dev/null | head
