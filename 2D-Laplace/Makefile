CC=mpicc
CFLAGS=-Wall -g -O0
LDLIBS=-lm
ALL=Init Solve InitPar Combine ReadGridParamsTest SolveParRow SolvePar

.PHONY: all
all: ${ALL}

Init: Init.o Utils.o 

Solve: Solve.o SerialSolver.o Utils.o  

InitPar: InitPar.o Utils.o

Combine: Combine.o Utils.o

Utils.o: Utils.h Utils.c

ReadGridParamsTest: ReadGridParamsTest.o Utils.o

SolveParRow: SolveParRow.o Utils.o

SolvePar: SolvePar.o Utils.o

.PHONY: clean
clean:
	rm -f *.o ${ALL}
