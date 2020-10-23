.PHONY: all clean run
CC=gcc
CXX=g++
CFLAGS=-std=c99 -O3 -march=native -ffast-math -funsafe-math-optimizations
LDFLAGS=-lm -lpthread -fopenmp
TARGET=train test
TRAIN_ITER=400

all: $(TARGET)

train: src/train.cpp src/lib/hmm_utils.cpp
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS) -Iinc

test: src/test.cpp src/lib/hmm_utils.cpp
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
