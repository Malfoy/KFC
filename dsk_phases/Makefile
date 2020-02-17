# example GATB-core makefile

.PHONY: all clean

# change this to the folder where gatb-core is
GATB=$(shell pwd)/gatb-core-stripped

GATB_LIB=$(GATB)/build/lib
CONFIG=$(GATB)/build/include/gatb/system/api/config.hpp

CXXFLAGS =  -std=c++0x  -I$(GATB)/src  -I$(GATB)/build/include -I$(GATB)/thirdparty
LDFLAGS=   -L$(GATB_LIB) -lgatbcore  -lpthread -lhdf5 -lz -std=c++0x -ldl 

ifeq ($(deb),1)
 CXXFLAGS+=-O0 -g
 else
 CXXFLAGS+=-O3
endif

#ifeq ($(gcov),1)
# CXXFLAGS+= -fprofile-arcs -ftest-coverage
#endif

SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,%,$(SRCS))

all: $(CONFIG) $(OBJS) 

$(CONFIG):
	cd gatb-core-stripped ; mkdir build ; cd build ; cmake .. ; make -j4

clean:
	rm -fv $(OBJS)

%: %.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LDFLAGS)

