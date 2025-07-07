SHELL = /bin/bash

CC ?= gcc
override CFLAGS += -std=gnu17 -O3 -march=native -I include -fopenmp -DNDEBUG

OBJ_SHARED = graph.o local_search.o chils_internal.o

OBJ_CHILS = main.o $(OBJ_SHARED)
OBJ_CHILS := $(addprefix bin/, $(OBJ_CHILS))

OBJ_LIB = chils.o $(OBJ_SHARED)
OBJ_LIB := $(addprefix bin/, $(OBJ_LIB))

DEP = $(OBJ_CHILS) $(OBJ_LIB)
DEP := $(sort $(DEP))

vpath %.c src
vpath %.h include

all : CHILS libCHILS.a

-include $(DEP:.o=.d)

CHILS : $(OBJ_CHILS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

libCHILS.a : $(OBJ_LIB)
	ar -rc $@ $^

bin/%.o : %.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@ $(LDFLAGS) $(LDLIBS)

.PHONY : clean
clean :
	rm -f CHILS libCHILS.a $(DEP) $(DEP:.o=.d)
