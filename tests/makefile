CC=g++
CFLAGS= -Wall -g -std=c++11 -I.. -O0
LDFLAGS=-flto -lpthread -g -fopenmp -lz  -Isparsepp  -flto ../smhasher/src/libSMHasherSupport.a

CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o) ../kfc.a
EXEC=tests
all: $(EXEC)



$(EXEC):   $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

../kfc.a: force
	make -C .. kfc.a

clean:
	rm -f *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)

force:
