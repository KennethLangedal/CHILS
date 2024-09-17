SHELL = /bin/bash

CC = gcc
CFLAGS = -std=gnu17 -O3 -march=native -I include -fopenmp -DNDEBUG

OBJ_PILS = main_pils.o graph.o reductions.o local_search.o pils.o tiny_solver.o \
		   neighborhood.o unconfined.o domination.o clique.o
OBJ_KERNEL = main_kernel.o graph.o reductions.o tiny_solver.o \
			 clique.o unconfined.o neighborhood.o domination.o \
			 single_edge.o extended_single_edge.o twin.o extended_twin.o heavy_vertex.o

OBJ_PILS := $(addprefix bin/, $(OBJ_PILS))
OBJ_KERNEL := $(addprefix bin/, $(OBJ_KERNEL))

DEP = $(OBJ_PILS) $(OBJ_KERNEL)
DEP := $(sort $(DEP))

vpath %.c src src/reductions
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
	rm -f PILS KERNEL $(DEP) $(DEP:.o=.d)
