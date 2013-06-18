#Compiler stuff
CC := gcc
FLAGS := -lm -g -pthread -fopenmp

SRC := clbm.c main.c gnuplot_i.c plot.c micro_bc.c fsi.c lattice.c macro_bc.c input.c output.c workerpool.c lyapunov.c flow.c
OBJ := $(SRC:.c=.o)
PROG := lbm

all: $(PROG) $(OBJ)

$(PROG): $(OBJ)
	$(CC) -o $@ $(OBJ) $(FLAGS)

%.o: %.c
	$(CC) -c $< -o $@ $(FLAGS)
	
clean:
	rm -rf $(OBJ) $(PROG)