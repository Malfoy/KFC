CC=g++
AR=gcc-ar

DEBUG ?= 1
ifeq ($(DEBUG), 1)
        CFLAGS+=-O0
        WARNS= -Wextra -Wpedantic -Wno-format -Wswitch-default -Wswitch-enum -Wfloat-equal -Wconversion -Wsign-conversion \
        -Wold-style-cast -Wuseless-cast -Wlogical-op -Wcast-align -Wtrampolines -Werror=enum-compare -Wstrict-aliasing=2 \
        -Werror=parentheses -Wnull-dereference -Werror=restrict -Werror=logical-op -Wsync-nand -Werror=main -Wshift-overflow=2 \
        -Wcast-qual -Werror=array-bounds -Werror=char-subscripts -Wshadow -Werror=ignored-qualifiers \
        -Werror=sequence-point -Werror=address -Wduplicated-branches -Wsign-compare -Wodr -Wnarrowing -Wsuggest-final-methods \
        -Wformat-signedness -Wrestrict -Werror=aggressive-loop-optimizations -Werror=missing-braces -Werror=uninitialized \
        -Wframe-larger-than=32768 -Werror=nonnull -Wno-unused-function -Werror=init-self -Werror=empty-body -Wdouble-promotion \
        -Wduplicated-cond -Werror=write-strings -Werror=return-type -Wredundant-decls \
        -Werror=volatile-register-var -Wsuggest-final-types
        DEBUG_SYMS=1
else
        CFLAGS+=-O3 -flto -march=native -mtune=native
        LDFLAGS+=-O3 -flto -march=native -mtune=native
        WARNS=-Wfatal-errors
endif

ASSERTS ?= $(DEBUG)
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
CFLAGS+=-std=c++17 -pipe -fopenmp ${WARNS}
LDFLAGS+=-lpthread -fopenmp -lz

INCS=-Ithirdparty/gatb-lite/include/ -isystem thirdparty/sparsepp -isystem thirdparty/BBHash -Ithirdparty/smhasher/src/ -Ithirdparty/robin-hood-hashing/src/include
EXT_BUILT_LIBS=thirdparty/smhasher/src/libSMHasherSupport.a
SUBMODULE_TOKEN=thirdparty/smhasher/README.md

CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o)
DEPS = $(OBJS:%.o=%.d)
KFC_OBJ = SolidSampler.o BitSet.o BloomFilter.o Hash.o CascadingBloomFilter.o index_min.o BitSet.o

EXEC=kfc_blue kfc_red kmerCountEvaluator
LIB=kfc_blue.a

all: $(EXEC) $(LIB) tests

kmerCountEvaluator: evaluator.o $(EXT_BUILT_LIBS)
	@echo "[LD] $@"
	@$(CC) -o $@ $^ $(LDFLAGS) $(EXT_BUILT_LIBS)

kfc_blue: $(KFC_OBJ) $(EXT_BUILT_LIBS) kfc_blue.o
	@echo "[LD] $@"
	@$(CC) -o $@ $^ $(LDFLAGS) $(EXT_BUILT_LIBS)

kfc_red: $(KFC_OBJ) $(EXT_BUILT_LIBS) kfc_red.o
	@echo "[LD] $@"
	@$(CC) -o $@ $^ $(LDFLAGS) $(EXT_BUILT_LIBS)

kfc_blue.a: $(KFC_OBJ)
	@echo "[AR] $@"
	@$(AR) rcs kfc_blue.a $(KFC_OBJ)

-include $(DEPS)

%.o: %.cpp $(SUBMODULE_TOKEN)
	@echo "[CC] $<"
	@$(CC) $(CFLAGS) $(INCS) -MMD -o $@ -c $<

thirdparty/smhasher/src/libSMHasherSupport.a: $(SUBMODULE_TOKEN)
	cd thirdparty/smhasher/src/ && cmake .
	$(MAKE) -sC thirdparty/smhasher/src/ SMHasherSupport

$(SUBMODULE_TOKEN):
	git submodule update --init

tests: kfc_blue.a
	@$(MAKE) -s -C tests/
	@echo "[run tests]"
	@tests/tests

# TODO release: we remove the 'kfc' binary file for convenience when tests; remove that later
clean:
	@echo "[clean]"
	@rm -f $(EXEC) $(LIB) $(OBJS) $(DEPS) kfc
	@$(MAKE) -s -C tests/ clean

rebuild: clean
	@$(MAKE) -s all

.PHONY: all tests clean rebuild
