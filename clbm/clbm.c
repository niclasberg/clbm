#include "clbm.h"
#include "fsi.h"
#include "macro_bc.h"
#include "micro_bc.h"
#include "lattice.h"
#include <stdio.h>
#include <stdlib.h>

/*
 * Public methods
 */

void lbm_init_state(FlowState * f_state, LbmState * lbm_state)
{
	unsigned int i, j, k, idx;
	double ux0, uy0, rho0;

	// Allocate space for the distributions
	for(i = 0; i < Q; ++i) {
		lbm_state->f[i] = (double *) malloc(f_state->lx* f_state->ly * sizeof(double));
		lbm_state->f_next[i] = (double *) malloc(f_state->lx* f_state->ly * sizeof(double));
	}

	// Initialize the distributions
	for(i = 0; i < f_state->lx; ++i) {
		for(j = 0; j < f_state->ly; ++j) {
			idx = i*f_state->ly + j;

			ux0 = f_state->u[0][idx];
			uy0 = f_state->u[1][idx];
			rho0 = f_state->rho[idx];

			for(k = 0; k < Q; ++k) {
				lbm_state->f[k][idx] = feq(k, rho0, ux0, uy0);
			}
		}
	}
}

void lbm_destroy_state(LbmState * lbm_state)
{
	unsigned int i;
	for(i = 0; i < Q; ++i) {
		free(lbm_state->f[i]);
		free(lbm_state->f_next[i]);
	}
}

void lbm_lattice_info()
{
	unsigned int i;
	printf("D2Q9 lattice\n");
	printf("i\tw\tcx\tcy\n");

	for(i = 0; i < Q; ++i) {
		printf("%d\t%.2f\t%d\t%d\n", (int)i, (double)weight[i], (int)cx[i], (int)cy[i]);
	}
}

void lbm_print_info()
{
	printf("Lbm parameters:\n");
}

void lbm_run(FlowState * f_state, LbmState * lbm_state)
{
	collide(f_state, lbm_state);
	stream(f_state, lbm_state);
	hydrovar(f_state, lbm_state);
	implement_bcs(f_state, lbm_state);
}

/*
 * Lbm implementation methods
 */
void collide(FlowState * f_state, LbmState * lbm_state)
{
	unsigned int i, j, k, idx;
	double omega = 1.0 / f_state->tau, gx0, gy0, ux0, uy0, rho0;

	for(i = 0; i < f_state->lx; ++i) {
		for(j = 0; j < f_state->ly; ++j) {
			idx = i*f_state->ly + j;

			gx0 = f_state->force[0][idx];
			gy0 = f_state->force[1][idx];
			rho0 = f_state->rho[idx];
			ux0 = f_state->u[0][idx];
			uy0 = f_state->u[1][idx];

			for(k = 0; k < Q; ++k)
				lbm_state->f[k][idx] = lbm_state->f[k][idx] -
						omega * (lbm_state->f[k][idx] - feq(k, rho0, ux0, uy0)) +
						3.0 * weight[k] * (cx[k]*gx0 + cy[k]*gy0);
		}
	}
}

void stream(FlowState * f_state, LbmState * lbm_state)
{
	unsigned int i, j, lx_m, lx_p, ly_m, ly_p, idx;
	for(i = 0; i < f_state->lx; ++i) {
		for(j = 0; j < f_state->ly; ++j) {
			idx = i*f_state->ly + j;

			// Permute the indices, corresponding to a periodic boundary condition
			// This also avoids segfaults :)
			lx_m = i==0 				? f_state->lx-1 : i-1;
			lx_p = i==(f_state->lx-1) 	? 0 			: i+1;
			ly_m = j==0 				? f_state->ly-1	: j-1;
			ly_p = j==(f_state->ly-1)	? 0				: j+1;

			lbm_state->f_next[0][idx] = lbm_state->f[0][idx];
			lbm_state->f_next[1][idx] = lbm_state->f[1][lx_m*f_state->ly + j];
			lbm_state->f_next[2][idx] = lbm_state->f[2][i*f_state->ly + ly_m];
			lbm_state->f_next[3][idx] = lbm_state->f[3][lx_p*f_state->ly + j];
			lbm_state->f_next[4][idx] = lbm_state->f[4][i*f_state->ly + ly_p];
			lbm_state->f_next[5][idx] = lbm_state->f[5][lx_m*f_state->ly + ly_m];
			lbm_state->f_next[6][idx] = lbm_state->f[6][lx_p*f_state->ly + ly_m];
			lbm_state->f_next[7][idx] = lbm_state->f[7][lx_p*f_state->ly + ly_p];
			lbm_state->f_next[8][idx] = lbm_state->f[8][lx_m*f_state->ly + ly_p];
		}
	}
}

void create_node(Node * node, unsigned int i, unsigned int j, FlowState * f_state, LbmState * lbm_state)
{
	unsigned int k, idx;

	idx = i*f_state->ly + j;

	// Copy coordinates
	node->coord[0] = i;
	node->coord[1] = j;

	// Copy density
	node->rho = f_state->rho[idx];

	// Copy force and velocity
	for(k = 0; k < DIM; ++k) {
		node->u[k] = f_state->u[k][idx];
		node->force[k] = f_state->force[k][idx];
	}

	// Copy the post-stream distribution functions
	for(k = 0; k < Q; ++k)
		node->f[k] = lbm_state->f_next[k][idx];
}

void hydrovar(FlowState * f_state, LbmState * lbm_state)
{
	unsigned int i, j;
	//unsigned int corner_count = 0;
	//Node * node = (Node *) malloc(sizeof(Node));
	//Node * corners[4];

	Node node;

	// Evaluate the hydrodynamic variables for all nodes
	// except the corners (they will usually need extrapolation
	// that depends on the value of the other nodes)
	for(i = 0; i < f_state->lx; ++i) {
		for(j = 0; j < f_state->ly; ++j) {
			create_node(&node, i, j, f_state, lbm_state);

			/*if(f_state->is_corner[i*f_state->ly + j] != 0) {
				// Save for later treatment
				corners[corner_count] = node;
				// And create a new node object to operate on
				node = (Node *) malloc(sizeof(Node));
				++corner_count;
			} else {*/
				eval_hydrovar(f_state, lbm_state, &node);
			//}
		}
	}

	/*free(node);

	// Evaluate for the corner nodes
	for(i = 0; i < corner_count; ++i) {
		eval_hydrovar(f_state, lbm_state, corners[i]);
		free(corners[i]);
	}*/
}

void eval_hydrovar(FlowState * f_state, LbmState * lbm_state, Node * node)
{
	unsigned int idx, k;
	idx = node->coord[0]*f_state->ly + node->coord[1];

	if(f_state->macro_bc[idx] != 0) {
		macro_bc(node, f_state, f_state->macro_bc[idx]);
	} else {
		node->u[0] = 0.0;
		node->u[1] = 0.0;
		node->rho = 0.0;

		for(k = 0; k < Q; ++k) {
			node->u[0] += cx[k] * node->f[k];
			node->u[1] += cy[k] * node->f[k];
			node->rho += node->f[k];
		}

		node->u[0] /= node->rho;
		node->u[1] /= node->rho;
	}

	// Copy back the state to the global arrays
	f_state->rho[idx] = node->rho;
	f_state->u[0][idx] = node->u[0];
	f_state->u[1][idx] = node->u[1];
}

void implement_bcs(FlowState * f_state, LbmState * lbm_state)
{
	unsigned int i, j, k, idx;
	Node node;

	for(i = 0; i < f_state->lx; ++i) {
		for(j = 0; j < f_state->ly; ++j) {
			idx = i*f_state->ly + j;

			if(f_state->micro_bc[idx] != 0) {
				create_node(&node, i, j, f_state, lbm_state);
				micro_bc(&node, f_state->micro_bc[idx]);

				for(k = 0; k < Q; ++k)
					lbm_state->f_next[k][idx] = node.f[k];
			}
		}
	}
}
