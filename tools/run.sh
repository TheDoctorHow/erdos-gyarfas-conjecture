#!/bin/sh
# run.sh <n> <girth> <shards> <tag>
N=$1; G=$2; S=$3; TAG=$4
mkdir -p out
i=0
while [ $i -lt $S ]; do
  EG_NEARFILE=out/near_${TAG}_$i.txt EG_TAG=${TAG}_$i \
    ./eghunt.exe $N $G o g m $i $S > out/log_${TAG}_$i.txt 2> out/err_${TAG}_$i.txt &
  i=$((i+1))
done
wait
echo "ALL SHARDS DONE $TAG"
cat out/near_${TAG}_*.txt 2>/dev/null | sort -u > out/near_${TAG}_ALL.txt
echo "near-misses: $(wc -l < out/near_${TAG}_ALL.txt)"
grep -h EGSTAT out/err_${TAG}_*.txt | awk '/seen=/{split($3,a,"=");split($4,b,"=");split($5,c,"=");s+=a[2];nn+=b[2];cc+=c[2]} END{printf "TOTAL seen=%d noC4C8=%d counterexamples=%d\n",s,nn,cc}'
grep -h COUNTEREXAMPLE out/log_${TAG}_*.txt 2>/dev/null | head
