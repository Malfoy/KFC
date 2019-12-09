# example GATB-core makefile

.PHONY: all clean

# change this to the folder where gatb-core is
GATB=$(shell pwd)/gatb-core-stripped

GATB_LIB=$(GATB)/build/lib
CONFIG=$(GATB)/build/src/gatb/system/api/config.hpp

CXXFLAGS =  -std=c++0x -O3  -I$(GATB)/src  -I$(GATB)/build/include -I$(GATB)/thirdparty 
LDFLAGS=   -L$(GATB_LIB) -lgatbcore  -lpthread -lhdf5 -lz -std=c++0x -ldl -static


SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,%,$(SRCS))

all: $(CONFIG) $(OBJS) 

$(CONFIG):
	cd gatb-core-stripped ; mkdir build ; cd build ; cmake .. ; make -j4

clean:
	rm -fv $(OBJS)

%: %.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LDFLAGS)

