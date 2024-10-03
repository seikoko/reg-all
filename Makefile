CXX = g++
CC = g++ -x c
LD = g++

DEBUG ?= 1
SAN ?= 0
COMMON_FLAGS := -I include -Wall -Wextra -Wconversion -ggdb3 -Wno-missing-field-initializers
CXXFLAGS = $(COMMON_FLAGS) -std=c++20
CFLAGS = $(COMMON_FLAGS)
LDFLAGS =

ifeq ($(DEBUG),1)
	CXXFLAGS += -Og
	CFLAGS += -Og
	LDFLAGS +=
else
	CXXFLAGS += -O2
	CFLAGS += -O2
	LDFLAGS +=
endif

ifeq ($(SAN),1)
	CXXFLAGS += -fsanitize=address,undefined
	CFLAGS   += -fsanitize=address,undefined
	LDFLAGS  += -fsanitize=address,undefined
endif

CXXFILES = $(shell find . -name '*.cxx')
CFILES   = $(shell find . -name '*.c')
CXXOBJ = $(patsubst ./src/%.cxx,bin/%.cxx.o,$(CXXFILES))
COBJ   = $(patsubst ./src/%.c,bin/%.c.o,$(CFILES))

bin/%.cxx.o: src/%.cxx
	$(CXX) -c -o $@ $< $(CXXFLAGS)

bin/%.c.o: src/%.c
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(CXXOBJ) $(COBJ)
	$(LD) -o $@ $^ $(LDFLAGS)

.PHONY: clean run
clean:
	rm -f main
	rm -f bin/*

run: main
	./$<

