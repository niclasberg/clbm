#include "macro_bc.h"
#include "data_types.h"
#include "macros.h"
#include <stdio.h>

int bc_east = 1, bc_north = 2, bc_west = 3, bc_south = 4,
	bc_north_east = 5, bc_north_west = 6, bc_south_west = 7, bc_south_east = 8;

void macro_bc(Node * node, FlowState * f_state, unsigned int type)
{
	if(type == bc_east)
		macro_bc_east(node, f_state);
	else if(type == bc_north)
		macro_bc_north(node, f_state);
	else if(type == bc_west)
		macro_bc_west(node, f_state);
	else if(type == bc_south)
		macro_bc_south(node, f_state);
	else if(type == bc_north_east)
		macro_bc_north_east(node, f_state);
	else if(type == bc_north_west)
		macro_bc_north_west(node, f_state);
	else if(type == bc_south_west)
		macro_bc_south_west(node, f_state);
	else if(type == bc_south_east)
		macro_bc_south_east(node, f_state);
	else
		printf("Unknown bc type: %d\n", type);
}

/*
 * Straight wall bcs:
 * Sets the velocity and computes the density from the known distribution functions
 */
void macro_bc_east(Node * node, FlowState * f_state)
{
	node->u[0] = f_state->u_ref * (-1 + 2 * (double)node->coord[1] / (f_state->ly - 1.0f));
	node->u[1] = 0.0;
	node->rho = (node->f[0] + 2*node->f[1] + node->f[2] + node->f[4] +
			2*node->f[5] + 2*node->f[8])/(node->u[0] + 1);
}

void macro_bc_west(Node * node, FlowState * f_state)
{
	node->u[0] = f_state->u_ref * (-1 + 2 *(double)node->coord[1] / (f_state->ly - 1.0f));
	node->u[1] = 0.0;
	node->rho = (node->f[0] + node->f[2] + 2*node->f[3] + node->f[4] +
			2*node->f[6] + 2*node->f[7])/(1 - node->u[0]);
}

void macro_bc_north(Node * node, FlowState * f_state)
{
	node->u[0] = f_state->u_ref;
	node->u[1] = 0.0;
	node->rho = (node->f[0] + node->f[1] + 2*node->f[2] + node->f[3] +
			2*node->f[5] + 2*node->f[6])/(node->u[1] + 1);
}

void macro_bc_south(Node * node, FlowState * f_state)
{
	node->u[0] = -f_state->u_ref;
	node->u[1] = 0.0;
	node->rho = (node->f[0] + node->f[1] + node->f[3] +
			2*node->f[4] + 2*node->f[7] + 2*node->f[8])/(1 - node->u[1]);
}

/*
 * Corner bcs:
 * Sets the velocity and extrapolates the density
 */
void macro_bc_north_east(Node * node, FlowState * f_state)
{
	node->u[0] = f_state->u_ref;
	node->u[1] = 0.0;
	node->rho = f_state->rho[(node->coord[0]-1)*f_state->ly + node->coord[1]  ] +
				f_state->rho[(node->coord[0]  )*f_state->ly + node->coord[1]-1] -
				f_state->rho[(node->coord[0]-1)*f_state->ly + node->coord[1]-1];
}

void macro_bc_north_west(Node * node, FlowState * f_state)
{
	node->u[0] = f_state->u_ref;
	node->u[1] = 0.0;
	node->rho = f_state->rho[(node->coord[0]+1)*f_state->ly + node->coord[1]  ] +
				f_state->rho[(node->coord[0]  )*f_state->ly + node->coord[1]-1] -
				f_state->rho[(node->coord[0]+1)*f_state->ly + node->coord[1]-1];
}

void macro_bc_south_west(Node * node, FlowState * f_state)
{
	node->u[0] = -f_state->u_ref;
	node->u[1] = 0.0;
	node->rho = f_state->rho[(node->coord[0]+1)*f_state->ly + node->coord[1]  ] +
				f_state->rho[(node->coord[0]  )*f_state->ly + node->coord[1]+1] -
				f_state->rho[(node->coord[0]+1)*f_state->ly + node->coord[1]+1];
}

void macro_bc_south_east(Node * node, FlowState * f_state)
{
	node->u[0] = -f_state->u_ref;
	node->u[1] = 0.0;
	node->rho = f_state->rho[(node->coord[0]-1)*f_state->ly + node->coord[1]  ] +
				f_state->rho[(node->coord[0]  )*f_state->ly + node->coord[1]+1] -
				f_state->rho[(node->coord[0]-1)*f_state->ly + node->coord[1]+1];
}
