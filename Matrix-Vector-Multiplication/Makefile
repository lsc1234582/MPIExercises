CC=mpicc
C_FLAGS=-Wall
LD_FLAGS=-lm

.PHONY: all
all: main serial

main: main.c
	$(CC) -o $@ $^ $(C_FLAGS) $(LD_FLAGS)

serial: serial.c
	$(CC) -o $@ $^ $(C_FLAGS) $(LD_FLAGS)

.PHONY: clean
clean:
	rm -f *.o main serial
