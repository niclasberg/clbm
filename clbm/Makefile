#Compiler stuff
CC := gcc
FLAGS := -lm -g -pthread -fopenmp -lgomp -O3

SRC := clbm.c main.c micro_bc.c fsi.c lattice.c macro_bc.c input.c output.c workerpool.c flow.c iohelpers.c
OBJ := $(SRC:.c=.o)
PROG := lbm
MEXPROG := lbm_ebf.mexa64
MEXSRC := lbm_ebf_matlab.c clbm.c micro_bc.c fsi.c lattice.c macro_bc.c input.c output.c flow.c iohelpers.c
BENCHSRC := clbm.c benchmark.c micro_bc.c fsi.c lattice.c macro_bc.c input.c output.c workerpool.c flow.c iohelpers.c
BENCHOBJ := $(BENCHSRC:.c=.o)
BENCHPROG := bench_prog

all: $(PROG) $(OBJ)

matlab: $(MEXPROG)

benchmark : $(BENCHPROG)
	./$(BENCHPROG)

$(BENCHPROG) : $(BENCHOBJ)
	$(CC) -o $@ $(BENCHOBJ) $(FLAGS)

$(MEXPROG) : $(MEXSRC)
	mex CC="gcc" CC="gcc" CFLAGS="\$$CFLAGS -fopenmp" LDFLAGS="\$$LDFLAGS -fopenmp" $(CFLAGS) $(MEXSRC) -o $(MEXPROG)

$(PROG): $(OBJ)
	$(CC) -o $@ $(OBJ) $(FLAGS)

%.o: %.c
	$(CC) -c $< -o $@ $(FLAGS)

clean:
	rm -rf $(OBJ) $(PROG) $(MEXPROG) $(BENCHPROG)