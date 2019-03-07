CC=g++
CFLAGS= -Wall -g -std=c++11 -I.. -I../thirdparty/smhasher/src/  -I../thirdparty/sparsepp -I../thirdparty/BBHash -I../thirdparty/ -O0
LDFLAGS=-lpthread -g -lz ../thirdparty/smhasher/src/libSMHasherSupport.a

CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o)
DEPS = $(OBJS:%.o=%.d)
EXEC=tests
all: $(EXEC)



$(EXEC): $(OBJS) ../kfc.a
	@echo "[LD] $@"
	@$(CC) -o $@ $^ $(LDFLAGS)

-include $(DEPS)

%.o: %.cpp
	@echo "[CC] $<"
	@$(CC) $(CFLAGS) -MMD -o $@ -c $<

../kfc.a:
	@echo "[build kfc.a]"
	@$(MAKE) -s -C .. kfc.a

clean:
	@echo "[clean tests/]"
	@rm -f $(OBJS) $(DEPS)


rebuild: clean
	@$(MAKE) -s all

