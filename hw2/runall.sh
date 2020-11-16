#!/bin/bash
source 0-activate.sh
echo "============ Running Preprocessing ============"
bash 1-preprocess.sh
echo "============= Extracting Features ============="
bash 2-extract-feat.sh
echo "================== Training ==================="
bash 3-train.sh
echo "=================== Testing ==================="
bash 4-test.sh
echo "Done!"
