CC=mpicc
CFLAGS=-Wall -g -O0
LDLIBS=-lm
ALL=Init InitPar Combine Solve SolveParRow SolvePar ReadGridParamsTest

.PHONY: all
all: ${ALL}

Init: Init.o Utils.o 

InitPar: InitPar.o Utils.o

Combine: Combine.o Utils.o

Solve: Solve.o Utils.o  

SolveParRow: SolveParRow.o Utils.o

SolvePar: SolvePar.o Utils.o

Utils.o: Utils.h Utils.c

ReadGridParamsTest: ReadGridParamsTest.o Utils.o

.PHONY: clean
clean:
	rm -f *.o ${ALL}
