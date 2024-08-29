SHELL = /bin/bash

CC = gcc
CFLAGS = -g -std=gnu17 -O3 -march=native -fopenmp -I include -D _GNU_SOURCE

OBJ = main.o graph.o local_search.o mwis.o reductions.o

OBJ := $(addprefix bin/, $(OBJ))

DEP = $(OBJ)
DEP := $(sort $(DEP))

vpath %.c src
vpath %.h include

all : MWIS

-include $(DEP:.o=.d)

MWIS : $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

bin/%.o : %.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@

.PHONY : clean
clean :
	rm -f MWIS
	rm -f $(DEP)
	rm -f $(DEP:.o=.d)
