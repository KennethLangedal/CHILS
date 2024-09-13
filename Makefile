SHELL = /bin/bash

CC = gcc
CFLAGS = -std=gnu17 -O3 -march=native -I include -fopenmp -DNDEBUG

OBJ_PILS = main_pils.o graph.o reductions.o local_search.o pils.o
OBJ_KERNEL = main_kernel.o graph.o reductions.o clique.o unconfined.o neighborhood.o domination.o

OBJ_PILS := $(addprefix bin/, $(OBJ_PILS))
OBJ_KERNEL := $(addprefix bin/, $(OBJ_KERNEL))

DEP = $(OBJ_PILS) $(OBJ_KERNEL)
DEP := $(sort $(DEP))

vpath %.c src
vpath %.c src/reductions
vpath %.h include

all : PILS KERNEL

-include $(DEP:.o=.d)

PILS : $(OBJ_PILS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

KERNEL : $(OBJ_KERNEL)
	$(CC) $(CFLAGS) -o $@ $^ -lm

bin/%.o : %.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@

.PHONY : clean
clean :
	rm -f PILS
	rm -f KERNEL
	rm -f $(DEP)
	rm -f $(DEP:.o=.d)
