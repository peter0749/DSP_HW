.PHONY: all clean run
CC=gcc
CXX=g++
CFLAGS=-std=c99 -O3 -march=native --fast-math
LDFLAGS=-lm -lpthread
TARGET=train test
TRAIN_ITER=200

all: $(TARGET)

lib/hmm_utils.o: src/lib/hmm_utils.c inc/lib/hmm_utils.h
	if [ ! -d ./lib ]; then mkdir -p ./lib; fi
	$(CC) -o $@ $< -c $(CFLAGS) -Iinc

train: src/train.cpp lib/hmm_utils.o
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS) -Iinc

test: src/test.cpp lib/hmm_utils.o
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS) -Iinc

clean:
	rm -f $(TARGET)
	rm -rf ./lib
	rm -f model_0*.txt result*.txt

train_all: train
	chmod u+x ./train_all.sh
	./train_all.sh ${TRAIN_ITER}

eval_all: all
	chmod u+x ./eval_all.sh
	./eval_all.sh
