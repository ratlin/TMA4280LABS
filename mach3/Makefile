CC=mpicc
CFLAGS=-fopenmp

TARGET=mach3
TARGETSRC=main.c
TARGETOBJ=$(patsubst %.c, %.o, $(TARGETSRC))

SRC=mach3.c
OBJ=$(patsubst %.c, %.o, $(SRC))

%.o: %.c
	$(CC) $(CFLAGS) -c $<

all: $(TARGET)

$(TARGET): $(TARGETOBJ) $(OBJ)
	$(CC) -o $@ $^ -lm -lgomp

clean:
	rm -f *.o $(TARGET) utest vtest

utest: utest.o $(OBJ)
	$(CC) -o $@ $^ -lm -lgomp

run_utest: utest
	OMP_NUM_THREADS=4 mpirun -n 4 ./utest

vtest: vtest.o $(OBJ)
	$(CC) -o $@ $^ -lm -lgomp

run_vtest: vtest
	OMP_NUM_THREADS=4 mpirun -n 8 ./vtest
