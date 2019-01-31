CC=mpicc
CFLAGS=-Wall -g -O0
LDLIBS=-lm
ALL=Init Solve InitParRow Combine ReadGridParamsTest

.PHONY: all
all: ${ALL}

Init: Init.o Utils.o 

Solve: Solve.o SerialSolver.o Utils.o  

InitParRow: InitParRow.o Utils.o

Combine: Combine.o Utils.o

ReadGridParamsTest: ReadGridParamsTest.o Utils.o

.PHONY: clean
clean:
	rm -f *.o ${ALL}
