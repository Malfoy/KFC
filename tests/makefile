CC=g++
CFLAGS= -Wall -g -std=c++11 -I..  -O0
LDFLAGS=-lpthread -g -lz

CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o)
DEPS = $(OBJS:%.o=%.d)
EXEC=tests
all: $(EXEC)



$(EXEC): $(OBJS) ../kfc_blue.a
	@echo "[LD] $@"
	@$(CC) -o $@ $^ $(LDFLAGS)

-include $(DEPS)

%.o: %.cpp
	@echo "[CC] $<"
	@$(CC) $(CFLAGS) -MMD -o $@ -c $<

../kfc_blue.a:
	@echo "[build kfc_blue.a]"
	@$(MAKE) -s -C .. kfc_blue.a

clean:
	@echo "[clean tests/]"
	@rm -f $(OBJS) $(DEPS)


rebuild: clean
	@$(MAKE) -s all

