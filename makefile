#~ CC=/usr/bin/g++
CC=g++
CFLAGS= -Wall -g -std=c++11   -pipe -funit-at-a-time -fopenmp -lz -Isparsepp
LDFLAGS= -lpthread -g -fopenmp -lz  -Isparsepp  -flto smhasher/src/libSMHasherSupport.a

CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o)


EXEC=kfc kmerCountEvaluator
# LIB=$(EXEC).a
all: $(EXEC)

kmerCountEvaluator:   evaluator.o
	$(CC) -o $@ $^ $(LDFLAGS)


kfc: kfc.o SolidSampler.o BitSet.o BloomFilter.o Hash.o	CascadingBloomFilter.o index_min.o
	$(CC) -o $@ $^ $(LDFLAGS)


%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -f *.o *.a
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
