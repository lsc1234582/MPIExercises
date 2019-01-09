CC=mpicc
C_FLAGS=-Wall
LD_FLAGS=-lm

.PHONY: all
all: main

main: main.c
	$(CC) -o main main.c $(C_FLAGS) $(LD_FLAGS)

.PHONY: clean
clean:
	rm -f *.o main
