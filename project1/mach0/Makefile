CC=gcc

TARGET=mach0
TARGETSRC=main.c
# substitutes .c with .o on $(TARGETSRC)
TARGETOBJ=$(patsubst %.c, %.o, $(TARGETSRC))

SRC=mach0.c
# substitutes .c with .o on $(SRC)
OBJ=$(patsubst %.c, %.o, $(SRC))

# target is anything.o, prereq is anything.c
# $< is the name of the first prereq
%.o: %.c
	$(CC) $(CFLAGS) -c $<

all: $(TARGET)

$(TARGET): $(TARGETOBJ) $(OBJ)
	$(CC) -o $@ $^ -lm

clean:
	rm -f *.o $(TARGET)

# $@ is the name of the target
# $^ is the names of all prereqs
utest: utest.o $(OBJ)
	$(CC) -o $@ $^ -lm

run_utest: utest
	./utest

# $@ is the name of the target
# $^ is the names of all prereqs
vtest: vtest.o $(OBJ)
	$(CC) -o $@ $^ -lm

run_vtest: vtest
	./vtest