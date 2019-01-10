CC=mpicc
C_FLAGS=-Wall -g -O3
LD_FLAGS=-lm

.PHONY: all
all: main

main: main.c
	$(CC) -o $@ $^ $(C_FLAGS) $(LD_FLAGS)

.PHONY: clean
clean:
	rm -f *.o main
