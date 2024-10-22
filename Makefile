SHELL = /bin/bash

CC = gcc
CFLAGS = -std=gnu17 -O3 -march=native -I include -fopenmp -DNDEBUG

OBJ_CHILS = main.o graph.o local_search.o chils.o
OBJ_CHILS := $(addprefix bin/, $(OBJ_CHILS))

DEP = $(OBJ_CHILS)
DEP := $(sort $(DEP))

vpath %.c src
vpath %.h include

all : CHILS

-include $(DEP:.o=.d)

CHILS : $(OBJ_CHILS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

bin/%.o : %.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@

.PHONY : clean
clean :
	rm -f CHILS $(DEP) $(DEP:.o=.d)
