CC=mpicc

TARGET=zeta1
TARGETSRC=main.c
TARGETOBJ=$(patsubst %.c, %.o, $(TARGETSRC))

SRC=zeta1.c
OBJ=$(patsubst %.c, %.o, $(SRC))


%.o: %.c
	$(CC) $(CFLAGS) -c $<

all: $(TARGET)

$(TARGET): $(TARGETOBJ) $(OBJ)
	$(CC) -o $@ $^ -lm

clean:
	rm -f *.o $(TARGET) utest vtest

utest: utest.o $(OBJ)
	$(CC) -o $@ $^ -lm

run_utest: utest
	mpirun -n 4 ./utest

vtest: vtest.o $(OBJ)
	$(CC) -o $@ $^ -lm

run_vtest: vtest
	mpirun -n 8 ./vtest