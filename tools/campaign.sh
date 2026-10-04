#!/bin/sh
# run from the folder that holds sweep.sh, egtool.exe and sh/snarkhunter-2.0b
# (build per RESULTS.md section 1; nauty and snarkhunter are NOT in this repo)
cd "$(dirname "$0")"
for N in 28 30; do
  echo "############ starting n=$N at $(date) ############"
  s=$(date +%s)
  ./sweep.sh $N 8
  echo "n=$N elapsed $(( $(date +%s)-s ))s"
done
