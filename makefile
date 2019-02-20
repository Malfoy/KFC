#~ CC=/usr/bin/g++
CC=g++

DEBUG ?= 1
ifeq ($(DEBUG), 1)
	CFLAGS+=-O0
	WARNS= -Wextra -Wpedantic -Wno-format -Wpadded -Wswitch-default -Wswitch-enum -Wfloat-equal -Wconversion -Wsign-conversion -Wold-style-cast -Wuseless-cast -Wlogical-op -Wcast-align -Wtrampolines -Werror=enum-compare -Wstrict-aliasing=2 -Werror=parentheses -Wnull-dereference -Werror=restrict -Werror=logical-op -Wsync-nand -Werror=main -Wshift-overflow=2 -Werror=pointer-sign -Wcast-qual -Werror=array-bounds -Werror=char-subscripts -Wshadow -Werror=ignored-qualifiers -Werror=sequence-point -Werror=address -Wduplicated-branches -Wsign-compare -Wodr -Wnarrowing -Wsuggest-final-methods  -Wformat-signedness -Wrestrict -Werror=aggressive-loop-optimizations -Werror=missing-braces -Werror=uninitialized -Wframe-larger-than=32768 -Werror=nonnull -Wno-unused-function -Werror=init-self -Werror=empty-body -Wdouble-promotion -Wfatal-errors -Werror=old-style-declaration -Wduplicated-cond -Werror=write-strings -Werror=return-type -Werror=volatile-register-var -Wsuggest-final-types -Werror=missing-parameter-type -Werror=implicit-int
	DEBUG_SYMS=1
else
	CFLAGS+=-O3 -flto -march=native -mtune=native
	LDFLAGS+=-O3 -flto -march=native -mtune=native
	WARNS=-Wfatal-errors
endif

ASSERTS ?= DEBUG
ifeq ($(ASSERTS), 1)
	CFLAGS+=-DDEBUG
else
	CFLAGS+=-DNDEBUG
endif

DEBUG_SYMS ?= 1
ifeq ($(DEBUG_SYMS), 1)
	CFLAGS+=-g
	LDFLAGS+=-g
endif

WARNS+= -Wall
CFLAGS+=-std=c++11 -pipe -fopenmp ${WARNS}
LDFLAGS+=-lpthread -fopenmp -lz

INCS=-Ithirdparty/gatb-lite/include/ -Ithirdparty/sparsepp -Ithirdparty/BBHash -Ithirdparty/smhasher/src/
LIBS=thirdparty/smhasher/src/libSMHasherSupport.a

CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o)
KFC_OBJ = kfc.o SolidSampler.o BitSet.o BloomFilter.o Hash.o CascadingBloomFilter.o index_min.o


EXEC=kfc kmerCountEvaluator
LIB=kfc.a
all: $(EXEC) $(LIB)

kmerCountEvaluator:   evaluator.o $(LIBS)
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)


kfc: $(KFC_OBJ) $(LIBS)
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

kfc.a: $(KFC_OBJ)
	ar rcs kfc.a $(KFC_OBJ)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS) $(INCS)

thirdparty/smhasher/src/libSMHasherSupport.a: FORCE
	cmake -S thirdparty/smhasher/src/ -B thirdparty/smhasher/src/
	$(MAKE) -C thirdparty/smhasher/src/ SMHasherSupport

clean:
	rm -f *.o *.a
	rm -rf $(EXEC)

rebuild: clean $(EXEC)

FORCE: ;
