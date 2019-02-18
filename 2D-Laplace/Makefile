CC=mpicc
CFLAGS=-Wall -g -O0
LDLIBS=-lm
ALL=Init Solve InitParRow Combine ReadGridParamsTest SolveParRow

.PHONY: all
all: ${ALL}

Init: Init.o Utils.o 

Solve: Solve.o SerialSolver.o Utils.o  

InitParRow: InitParRow.o Utils.o

Combine: Combine.o Utils.o

Utils.o: Utils.h Utils.c

ReadGridParamsTest: ReadGridParamsTest.o Utils.o

SolveParRow: SolveParRow.o Utils.o

.PHONY: clean
clean:
	rm -f *.o ${ALL}
