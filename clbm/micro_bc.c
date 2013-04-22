#include "micro_bc.h"
#include "lattice.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int bc_zou_he_east = 1, bc_zou_he_north = 2, bc_zou_he_west = 3, bc_zou_he_south = 4,
	bc_zou_he_north_east = 5, bc_zou_he_north_west = 6, bc_zou_he_south_west = 7, bc_zou_he_south_east = 8,
	bc_regularized_east = 9, bc_regularized_north = 10, bc_regularized_west = 11, bc_regularized_south = 12,
	bc_diffusive_east = 13, bc_diffusive_north = 14, bc_diffusive_west = 15, bc_diffusive_south = 16;

int INIT_DIFFUSIVE = 1, INIT_ZOU_HE = 2, INIT_REGULARIZED = 4;

unsigned int nx, ny;
double * f_prev[Q];

void bc_init(int types, LbmState * lbm_state, unsigned int lx, unsigned int ly)
{
	unsigned int i, k, idx;
	if(types & INIT_DIFFUSIVE) {
		nx = lx;
		ny = ly;
		for(k = 0; k < Q; ++k) {
			f_prev[k] = (double *) malloc(lx*ly*sizeof(double));

			for(i = 0; i < lx*ly; ++i) {
				f_prev[k][i] = lbm_state->f[k][i];
			}
		}
	}
}

void bc_destroy()
{
	if(f_prev[0]) {
		unsigned int k;
		for(k = 0; k < Q; ++k) {
			free(f_prev[k]);
		}
	}
}

void micro_bc(Node * node, unsigned int type)
{
	if(type == bc_zou_he_east)
		zou_he_east(node);
	else if(type == bc_zou_he_north)
		zou_he_north(node);
	else if(type == bc_zou_he_west)
		zou_he_west(node);
	else if(type == bc_zou_he_south)
		zou_he_south(node);
	else if(type == bc_zou_he_north_east)
		zou_he_north_east(node);
	else if(type == bc_zou_he_north_west)
		zou_he_north_west(node);
	else if(type == bc_zou_he_south_west)
		zou_he_south_west(node);
	else if(type == bc_zou_he_south_east)
		zou_he_south_east(node);
	else if(type == bc_regularized_east)
		bc_regularized_straight(node->f, node->rho, node->u[0], node->u[1], 0);
	else if(type == bc_regularized_north)
		bc_regularized_straight(node->f, node->rho, node->u[0], node->u[1], 1);
	else if(type == bc_regularized_west)
		bc_regularized_straight(node->f, node->rho, node->u[0], node->u[1], 2);
	else if(type == bc_regularized_south)
		bc_regularized_straight(node->f, node->rho, node->u[0], node->u[1], 3);
	else if(type == bc_diffusive_north)
		diffusive_north(node);
	else if(type == bc_diffusive_south)
		diffusive_south(node);
	else
		printf("Unknown micro_bc type: %d\n", type);
}

void zou_he_east(Node * node)
{
	node->f[3] = node->f[1] - (2*node->rho*node->u[0])/3.0;
	node->f[6] = node->f[8] + 0.5*(node->f[4] - node->f[2]) - (node->rho*node->u[0])/6.0 + (node->rho*node->u[1])/2.0;
	node->f[7] = node->f[5] + 0.5*(node->f[2] - node->f[4]) - (node->rho*node->u[0])/6.0 - (node->rho*node->u[1])/2.0;
}

void zou_he_west(Node * node)
{
	node->f[1] = node->f[3] + (2.0*node->rho*node->u[0])/3.0;
	node->f[5] = node->f[7] + 0.5*(node->f[4] - node->f[2]) + (node->rho*node->u[0])/6.0 + (node->rho*node->u[1])/2.0;
	node->f[8] = node->f[6] + 0.5*(node->f[2] - node->f[4]) + (node->rho*node->u[0])/6.0 - (node->rho*node->u[1])/2.0;
}

void zou_he_north(Node * node)
{
	node->f[4] = node->f[2] - (2.0*node->rho*node->u[1])/3.0;
	node->f[7] = node->f[5] + 0.5*(node->f[1] - node->f[3]) - (node->rho*node->u[0])/2.0 - (node->rho*node->u[1])/6.0;
	node->f[8] = node->f[6] + 0.5*(node->f[3] - node->f[1]) + (node->rho*node->u[0])/2.0 - (node->rho*node->u[1])/6.0;
}

void zou_he_south(Node * node)
{
	node->f[2] = node->f[4] + (2*node->rho*node->u[1])/3;
	node->f[5] = node->f[7] + 0.5*(node->f[3] - node->f[1]) + (node->rho*node->u[0])/2.0 + (node->rho*node->u[1])/6.0;
	node->f[6] = node->f[8] + 0.5*(node->f[1] - node->f[3]) - (node->rho*node->u[0])/2.0 + (node->rho*node->u[1])/6.0;
}

/*
 * Corner boundary conditions
 */
void zou_he_north_east(Node * node)
{
	node->f[3] = node->f[1] - (2*node->rho*node->u[0])/3;
	node->f[4] = node->f[2] - (2*node->rho*node->u[1])/3;
	node->f[6] = node->rho/2 - node->f[1] - node->f[2] - node->f[5] - node->f[0]/2 +
				(node->rho*node->u[0])/3 + (node->rho*node->u[1])/2;
	node->f[7] = node->f[5] - (node->rho*node->u[0])/6 - (node->rho*node->u[1])/6;
	node->f[8] = node->rho/2 - node->f[1] - node->f[2] - node->f[5] - node->f[0]/2 +
				(node->rho*node->u[0])/2 + (node->rho*node->u[1])/3;
}

void zou_he_north_west(Node * node)
{
	node->f[1] = node->f[3] + (2*node->rho*node->u[0])/3;
	node->f[4] = node->f[2] - (2*node->rho*node->u[1])/3;
	node->f[5] = node->rho/2 - node->f[2] - node->f[3] - node->f[6] - node->f[0]/2 -
				(node->rho*node->u[0])/3 + (node->rho*node->u[1])/2;
	node->f[7] = node->rho/2 - node->f[2] - node->f[3] - node->f[6] - node->f[0]/2 -
				(node->rho*node->u[0])/2 + (node->rho*node->u[1])/3;
	node->f[8] = node->f[6] + (node->rho*node->u[0])/6 - (node->rho*node->u[1])/6;
}

void zou_he_south_east(Node * node)
{
	node->f[2] = node->f[4] + (2*node->rho*node->u[1])/3;
	node->f[3] = node->f[1] - (2*node->rho*node->u[0])/3;
	node->f[5] = node->rho/2 - node->f[1] - node->f[4] - node->f[8] - node->f[0]/2 +
				(node->rho*node->u[0])/2 - (node->rho*node->u[1])/3;
	node->f[6] = node->f[8] - (node->rho*node->u[0])/6 + (node->rho*node->u[1])/6;
	node->f[7] = node->rho/2 - node->f[1] - node->f[4] - node->f[8] - node->f[0]/2 +
				(node->rho*node->u[0])/3 - (node->rho*node->u[1])/2;
}

void zou_he_south_west(Node * node)
{
	node->f[1] = node->f[3] + (2*node->rho*node->u[0])/3;
	node->f[2] = node->f[4] + (2*node->rho*node->u[1])/3;
	node->f[5] = node->f[7] + (node->rho*node->u[0])/6 + (node->rho*node->u[1])/6;
	node->f[6] = node->rho/2 - node->f[3] - node->f[4] - node->f[7] - node->f[0]/2 -
				(node->rho*node->u[0])/2 - (node->rho*node->u[1])/3;
	node->f[8] = node->rho/2 - node->f[3] - node->f[4] - node->f[7] - node->f[0]/2 -
				(node->rho*node->u[0])/3 - (node->rho*node->u[1])/2;
}

/*
 * Regularized bcs:
 * 	Evaluates the off-equilibrium parts of all distributions, applies non-equilibrium
 * 	bounce back to find the unknown distributions and uses the result to
 * 	find the stress tensor at the boundary.
 * 	This value is then used to evaluate all the distributions.
 * 	See Jonas Latt: Boundary review for more info.
 */

// Helper function to permute indices
unsigned int perm4(unsigned int i, unsigned int j)
{
	if(i <= 4)
		return 1 + (i+j-1)%4;
	else
		return 5 + (i+j-5)%4;
}

void bc_regularized_straight(double * f, double rho, double u, double v, int dir)
{
	unsigned int k;
	double fneq[Q];
	double Pi_xx, Pi_yy, Pi_xy, Q_xx, Q_xy, Q_yy;

	// Compute non-equilibrium part of the distributions
	for(k=0; k < Q; ++k)
		fneq[k] = f[k] - feq(k, rho, u, v);

	// Bounce-back on unknown distributions
	fneq[perm4(3, dir)] = fneq[perm4(1, dir)];
	fneq[perm4(6, dir)] = fneq[perm4(8, dir)];
	fneq[perm4(7, dir)] = fneq[perm4(5, dir)];

	// Evaluate stress tensor
	Pi_xx = Pi_xy = Pi_yy = 0;
	for(k=0; k < Q; ++k) {
		Pi_xx += cx[k]*cx[k]*fneq[k];
		Pi_xy += cx[k]*cy[k]*fneq[k];
		Pi_yy += cy[k]*cy[k]*fneq[k];
	}

	// Compute resulting distribution functions
	for(k=0; k < Q; ++k) {
		Q_xx = cx[k]*cx[k] - cs2;
		Q_xy = cx[k]*cy[k];
		Q_yy = cy[k]*cy[k] - cs2;

		f[k] = feq(k, rho, u, v) +
			weight[k]/(2*cs2*cs2) * (Q_xx*Pi_xx + 2*Q_xy*Pi_xy + Q_yy*Pi_yy);
	}
}

/*
 * Diffusive Bcs
 */
void diffusive_south(Node * node)
{
	unsigned int idx = node->coord[0] * ny + node->coord[1];
	double L;

	node->f[4] = 0.5 * (node->f[4] + f_prev[4][idx]);
	node->f[7] = 0.5 * (node->f[7] + f_prev[7][idx]);
	node->f[8] = 0.5 * (node->f[8] + f_prev[8][idx]);

	L = (node->f[4] + node->f[7] + node->f[8]) /
		(feq(4, node->rho, node->u[0], node->u[1]) + feq(7, node->rho, node->u[0], node->u[1]) + feq(8, node->rho, node->u[0], node->u[1]));

	node->f[2] = L * feq(2, node->rho, node->u[0], node->u[1]);
	node->f[5] = L * feq(5, node->rho, node->u[0], node->u[1]);
	node->f[6] = L * feq(6, node->rho, node->u[0], node->u[1]);

	f_prev[4][idx] = node->f[4];
	f_prev[7][idx] = node->f[7];
	f_prev[8][idx] = node->f[8];
}

void diffusive_north(Node * node)
{
	unsigned int idx = node->coord[0] * ny + node->coord[1];
	double L;

	node->f[2] = 0.5 * (node->f[2] + f_prev[2][idx]);
	node->f[5] = 0.5 * (node->f[5] + f_prev[5][idx]);
	node->f[6] = 0.5 * (node->f[6] + f_prev[6][idx]);

	L = (node->f[2] + node->f[5] + node->f[6]) /
		(feq(2, node->rho, node->u[0], node->u[1]) + feq(5, node->rho, node->u[0], node->u[1]) + feq(6, node->rho, node->u[0], node->u[1]));

	node->f[4] = L * feq(4, node->rho, node->u[0], node->u[1]);
	node->f[7] = L * feq(7, node->rho, node->u[0], node->u[1]);
	node->f[8] = L * feq(8, node->rho, node->u[0], node->u[1]);

	f_prev[2][idx] = node->f[2];
	f_prev[5][idx] = node->f[5];
	f_prev[6][idx] = node->f[6];
}
