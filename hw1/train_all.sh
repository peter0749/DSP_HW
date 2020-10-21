#!/bin/bash
NUM_JOBS=8
for ((i=1; i<=5; ++i)); do
    echo "${1} model_init.txt data/train_seq_0${i}.txt model_0${i}.txt"
done | xargs -P${NUM_JOBS} -n4 ./train
