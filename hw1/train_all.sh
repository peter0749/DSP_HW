#!/bin/bash
for ((i=1; i<=5; ++i)); do
    ./train ${1} model_init.txt data/train_seq_0${i}.txt model_0${i}.txt
done
