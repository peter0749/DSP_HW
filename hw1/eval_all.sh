#!/bin/bash
./test modellist.txt data/test_seq.txt result.txt
awk '{ print $1 }' result.txt > result_label_only.txt
WRONG_N=$(diff data/test_lbl.txt result_label_only.txt | grep "^>" | wc -l)
TOTAL=$(wc -l data/test_lbl.txt | awk '{print $1}')
CORRECT=$((TOTAL - WRONG_N))
printf "%.2f%%\n" $(bc <<< "scale=4; $CORRECT / $TOTAL * 100")
