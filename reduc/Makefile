CC=mpicc
TARGET=recursive_doubling_sum

all:
	$(CC) -c recursive_doubling_sum.c
	$(CC) -c main.c
	
	$(CC) -o $(TARGET) recursive_doubling_sum.o main.o  -lm

clean:
	rm -f *.o $(TARGET)