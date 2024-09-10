SHELL = /bin/bash

CC = gcc
CFLAGS = -std=gnu17 -O3 -march=native -I include -fopenmp -DNDEBUG

OBJ = main.o graph.o reductions.o local_search.o pils.o

OBJ := $(addprefix bin/, $(OBJ))

DEP = $(OBJ)
DEP := $(sort $(DEP))

vpath %.c src
vpath %.h include

all : PILS

-include $(DEP:.o=.d)

PILS : $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

bin/%.o : %.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@

.PHONY : clean
clean :
	rm -f PILS
	rm -f $(DEP)
	rm -f $(DEP:.o=.d)
