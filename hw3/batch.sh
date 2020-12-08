#/bin/bash
set -x
source setup.sh

TIMELIMIT="60s"

make
make map

perl separator_big5.pl corpus.txt > corpus.seg
ngram-count -text corpus.seg -write count.2 -order 2
ngram-count -text corpus.seg -write count.3 -order 3
ngram-count -read count.2 -lm bigram.lm -order 2 -unk
ngram-count -read count.3 -lm trigram.lm -order 3 -unk

if [[ 1 -eq 1 ]]; then
    cat input_list.txt | while read line; do
        fname=${line%.*}
        perl separator_big5.pl $line > $fname.seg;
        timeout "$TIMELIMIT" ./mydisambig $fname.seg ZhuYin-Big5.map bigram.lm $fname.mybigram.out
        timeout "$TIMELIMIT" ./mydisambig $fname.seg ZhuYin-Big5.map trigram.lm $fname.mytrigram.out
        timeout "$TIMELIMIT" disambig -text $fname.seg -map ZhuYin-Big5.map -lm bigram.lm -order 2 > $fname.bigram.out
        timeout "$TIMELIMIT" disambig -text $fname.seg -map ZhuYin-Big5.map -lm trigram.lm -order 3 > $fname.trigram.out
    done
fi
