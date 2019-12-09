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

WARNS+= -Wall
CFLAGS+=-std=c++14 -pipe -fopenmp ${WARNS}
LDFLAGS+=-lpthread -fopenmp -lz

INCS=-isystem thirdparty/gatb-lite/include/ -isystem thirdparty/robin-hood-hashing/src/include
SUBMODULE_TOKEN=thirdparty/gatb-lite/CMakeLists.txt

CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o)
DEPS = $(OBJS:%.o=%.d)
KFC_BLUE_OBJ = SolidSampler.o BitSet.o BloomFilter.o Hash.o CascadingBloomFilter.o BitSet.o

EXEC=kfc_blue kfc_red kmerCountEvaluator
LIB=kfc_blue.a

all: $(EXEC) $(LIB) tests

kmerCountEvaluator: evaluator.o
	@echo "[LD] $@"
	@mkdir --parents ./bin
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)
	@mv kmerCountEvaluator bin/

kfc_blue: $(KFC_BLUE_OBJ) $(EXT_BUILT_LIBS) kfc_blue.o
	@echo "[LD] $@"
	@mkdir --parents ./bin
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)
	@mv kfc_blue bin/

kfc_red: kfc_red.o
	@echo "[LD] $@"
	@mkdir --parents ./bin
	@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)
	@mv kfc_red bin/

kfc_blue.a: $(KFC_BLUE_OBJ)
	@echo "[AR] $@"
	@$(AR) rcs kfc_blue.a $(KFC_BLUE_OBJ)

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

# TODO release: we remove the 'kfc' binary file for convenience when tests; remove that later
clean:
	@echo "[clean]"
	@rm -f $(EXEC) $(LIB) $(OBJS) $(DEPS) kfc
	@$(MAKE) -s -C tests/ clean

rebuild: clean
	@$(MAKE) -s all

check_buildsys: $(SUBMODULE_TOKEN)
	$(CC) --version
	$(AR) --version
	@echo CFLAGS=$(CFLAGS)
	@echo LDFLAGS=$(LDFLAGS)

.PHONY: all tests clean rebuild
