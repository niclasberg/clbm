#ifndef MACRO_BC_H_
#define MACRO_BC_H_
#include "data_types.h"

int bc_east, bc_north, bc_west, bc_south, bc_north_east, bc_north_west, bc_south_west, bc_south_east;
void macro_bc(Node *, FlowState *, unsigned int);

static void macro_bc_east(Node *, FlowState *);
static void macro_bc_west(Node *, FlowState *);
static void macro_bc_north(Node *, FlowState *);
static void macro_bc_south(Node *, FlowState *);
static void macro_bc_north_east(Node *, FlowState *);
static void macro_bc_north_west(Node *, FlowState *);
static void macro_bc_south_west(Node *, FlowState *);
static void macro_bc_south_east(Node *, FlowState *);


#endif /* MACRO_BC_H_ */
