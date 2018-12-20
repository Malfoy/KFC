#~ CC=/usr/bin/g++
CC=g++
CFLAGS= -Wall -g -std=c++11  -flto -pipe -funit-at-a-time -fopenmp -lz -Isparsepp -flto
LDFLAGS=-flto -lpthread -g -fopenmp -lz  -Isparsepp  -flto smhasher/src/libSMHasherSupport.a

CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o)


EXEC=kfc
all: $(EXEC)



kfc:   $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -f *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
