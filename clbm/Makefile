#Compiler stuff
CC := gcc
FLAGS := -lm -g -pthread

SRC := clbm.c main.c micro_bc.c fsi.c lattice.c macro_bc.c input.c output.c workerpool.c flow.c iohelpers.c
OBJ := $(SRC:.c=.o)
PROG := lbm
MEXPROG := lbm_ebf.mexa64
MEXSRC := lbm_ebf_matlab.c clbm.c micro_bc.c fsi.c lattice.c macro_bc.c input.c output.c flow.c iohelpers.c

all: $(PROG) $(OBJ)

matlab: $(MEXPROG)

$(MEXPROG) : $(MEXSRC)
	mex CC="gcc" $(CFLAGS) $(MEXSRC) -o $(MEXPROG)

$(PROG): $(OBJ)
	$(CC) -o $@ $(OBJ) $(FLAGS)

%.o: %.c
	$(CC) -c $< -o $@ $(FLAGS)
	
clean:
	rm -rf $(OBJ) $(PROG) $(MEXPROG)