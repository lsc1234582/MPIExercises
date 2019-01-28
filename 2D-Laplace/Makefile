CC=mpicc
CFLAGS=-Wall -g -O3
LDLIBS=-lm

.PHONY: all
all: InitCondition Solver

InitCondition: InitCondition.c Utils.h

Solver: Solver.c Utils.h

.PHONY: clean
clean:
	rm -f *.o InitCondition Solver
