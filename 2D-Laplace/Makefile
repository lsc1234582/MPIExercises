CC=mpicc
CFLAGS=-Wall -g -O3
LDLIBS=-lm
ALL=Init Solve

.PHONY: all
all: ${ALL}

Init: Init.o Utils.o 

Solve: Solve.o SerialSolver.o Utils.o  

.PHONY: clean
clean:
	rm -f *.o ${ALL}
