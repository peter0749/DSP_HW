#!/bin/bash
source 0-activate.sh
mkdir -p "./log"

echo "========= Training all combinations ========="
awk '{ print $1 " " $2 " " $3 " " $4 }' "$1" | sort -n | uniq | sort -R | xargs -n4 -P$2 ./3-train-search.sh 2>&1 > "./log/search.train.log"

echo "========== Testing all combinations ========="
sort -n "$1" | uniq | sort -R | xargs -n6 -P$2 ./4-test-search.sh 2>&1 > "./log/search.test.log"

echo "================= Summarizing ==============="
find viterbi/mono -type f | grep -P "viterbi/mono/param-.*/test.macc" | xargs cat | sort -r
