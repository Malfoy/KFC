
CC=g++
AR=gcc-ar

DEBUG ?= 0
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
        CFLAGS+=-O3 -fno-fat-lto-objects -flto=jobserver -march=native -mtune=native
        LDFLAGS+=-fuse-linker-plugin
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

ifeq ($(VERBOSE), 1)
	SHELL=sh -x
endif

WARNS+= -Wall
CFLAGS+=-std=c++14 -pipe -fopenmp -mavx2 ${WARNS}
LDFLAGS+=-lpthread -fopenmp -lz

INCS=-isystem thirdparty/gatb-lite/include/ -isystem thirdparty/robin-hood-hashing/src/include
SUBMODULE_TOKEN=thirdparty/gatb-lite/CMakeLists.txt

CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o)
DEPS = $(OBJS:%.o=%.d)
KFC_BLUE_OBJ = SolidSampler.o BitSet.o BloomFilter.o Hash.o CascadingBloomFilter.o BitSet.o
KFC_RED_OBJ = 

EXEC=kfc_blue kfc_red kmerCountEvaluator SuperKmerCount
LIB=kfc_blue.a

all: $(EXEC) $(LIB) tests

bin:
	@mkdir --parent ./bin

kmerCountEvaluator: bin bin/kmerCountEvaluator
bin/kmerCountEvaluator: evaluator.o
	@echo "[LD] $@"
	+@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

kfc_blue: bin bin/kfc_blue
bin/kfc_blue: $(KFC_BLUE_OBJ) kfc_blue.o
	@echo "[LD] $@"
	+@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

kfc_red: bin bin/kfc_red
bin/kfc_red: $(KFC_RED_OBJ) kfc_red.o
	@echo "[LD] $@"
	+@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

SKC_tests: bin bin/SKC_tests
bin/SKC_tests: $(KFC_RED_OBJ) SKC_tests.o
	@echo "[LD] $@"
	+@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

kfc_blue.a: $(KFC_BLUE_OBJ)
	@echo "[AR] $@"
	@$(AR) rcs $@ $(KFC_BLUE_OBJ)

-include $(DEPS)

%.o: %.cpp $(SUBMODULE_TOKEN)
	@echo "[CC] $<"
	@$(CC) $(CFLAGS) $(INCS) -MMD -o $@ -c $<

$(SUBMODULE_TOKEN):
	git submodule update --init

tests: kfc_blue.a
	@$(MAKE) -s -C tests/
	@echo "[run tests]"
	@tests/tests

# TODO release: we remove the 'kfc' and 'bin/kfc' binary files for convenience when tests; remove that later
clean:
	@echo "[clean]"
	@rm -rf bin $(LIB) $(OBJS) $(DEPS) kfc
	@$(MAKE) -s -C tests/ clean

rebuild: clean
	@$(MAKE) -s all

check_buildsys: $(SUBMODULE_TOKEN)
	$(CC) --version
	$(AR) --version
	@echo CFLAGS=$(CFLAGS)
	@echo LDFLAGS=$(LDFLAGS)

.PHONY: all $(EXEC) tests clean rebuild
