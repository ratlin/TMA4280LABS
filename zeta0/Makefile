CC=gcc

TARGET=zeta0
TARGETSRC=main.c
TARGETOBJ=$(patsubst %.c, %.o, $(TARGETSRC))

SRC=zeta0.c
OBJ=$(patsubst %.c, %.o, $(SRC))


%.o: %.c
	$(CC) $(CFLAGS) -c $<

all: $(TARGET)

$(TARGET): $(TARGETOBJ) $(OBJ)
	$(CC) -o $@ $^ -lm

clean:
	rm -f *.o $(TARGET)

utest: utest.o $(OBJ)
	$(CC) -o $@ $^ -lm

run_utest: utest
	./utest

vtest: vtest.o $(OBJ)
	$(CC) -o $@ $^ -lm

run_vtest: vtest
	./vtest