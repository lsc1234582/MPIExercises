CC=mpicc

ifeq (${DEBUG}, 0)
	CFLAGS=-Wall -g -O3 -DNDEBUG
else
	CFLAGS=-Wall -g -O0 -DDEBUG
endif
LDLIBS=-lm
ALL=Init InitPar Combine Solve SolveParRow SolvePar ReadGridParamsTest Compare

.PHONY: all
all: ${ALL}

Init: Init.o Utils.o

InitPar: InitPar.o Utils.o

Combine: Combine.o Utils.o

Solve: Solve.o Utils.o

SolveParRow: SolveParRow.o Utils.o

SolvePar: SolvePar.o Utils.o

ReadGridParamsTest: ReadGridParamsTest.o Utils.o

Compare: Compare.o Utils.o

Utils.o: Utils.h Utils.c

.PHONY: clean
clean:
	rm -f *.o ${ALL}
