#include "clbm.h"
#include "fsi.h"
#include "macro_bc.h"
#include "micro_bc.h"
#include "lattice.h"
#include <stdio.h>
#include <stdlib.h>
#include "iohelpers.h"
#include <math.h>

/* Forward declarations */
void collide(FlowState * , LbmState *);
void stream(FlowState * , LbmState *);
void hydrovar_and_bcs(FlowState * , LbmState *);
void swap_states(LbmState *);

/*
 * Public methods
 */

LbmState * lbm_alloc_state(unsigned int lx, unsigned int ly)
{
	unsigned int i;
	LbmState * lbm_state = malloc(sizeof(LbmState));

	/* Domain size */
	lbm_state->lx = lx;
	lbm_state->ly = ly;

	/* Allocate space for the distributions */
	for(i = 0; i < Q; ++i) {
		lbm_state->f[i] = (double *) malloc(lx * ly * sizeof(double));
		lbm_state->f_next[i] = (double *) malloc(lx * ly * sizeof(double));
	}

	return lbm_state;
}

void lbm_free_state(LbmState * lbm_state)
{
	unsigned int i;
	for(i = 0; i < Q; ++i) {
		free(lbm_state->f[i]);
		lbm_state->f[i] = NULL;
		free(lbm_state->f_next[i]);
		lbm_state->f_next[i] = NULL;
	}
	lbm_state->lx = 0;
	lbm_state->ly = 0;

	free(lbm_state);
}

void lbm_init_state(FlowState * f_state, LbmState * lbm_state)
{
	unsigned int i, j, k, idx;
	double ux0, uy0, rho0;

	/* Initialize the distributions */
	for(i = 0; i < f_state->lx; ++i) {
		for(j = 0; j < f_state->ly; ++j) {
			idx = i*f_state->ly + j;

			ux0 = f_state->u[0][idx];
			uy0 = f_state->u[1][idx];
			rho0 = f_state->rho[idx];

			for(k = 0; k < Q; ++k) {
				lbm_state->f[k][idx] = feq(k, rho0, ux0, uy0);
				lbm_state->f_next[k][idx] = lbm_state->f[k][idx];
			}
		}
	}
}

LbmState * lbm_clone_state(const LbmState * src)
{
	LbmState * dest = lbm_alloc_state(src->lx, src->ly);
	lbm_copy_state(src, dest);
	return dest;
}

void lbm_copy_state(const LbmState * src, LbmState * dest)
{
	unsigned int i, j;
	#pragma omp for private(i)
	for(j = 0; j < dest->lx*dest->ly; ++j) {
		for(i = 0; i < Q; ++i) {
			dest->f[i][j] = src->f[i][j];
			dest->f_next[i][j] = src->f_next[i][j];
		}
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

void lbm_run(FlowState * f_state, LbmState * lbm_state)
{
	#pragma omp parallel
	{
		collide(f_state, lbm_state);
		stream(f_state, lbm_state);
		hydrovar_and_bcs(f_state, lbm_state);
	}
	swap_states(lbm_state);
}

/*
 * Lbm implementation methods
 */
void collide(FlowState * f_state, LbmState * lbm_state)
{
	unsigned int i, k;
	unsigned int nodes = lbm_state->lx * lbm_state->ly;
	double omega = 1.0 / f_state->tau, gx0, gy0, ux0, uy0, rho0;

	#pragma omp for private(i, k, gx0, gy0, ux0, uy0, rho0)
	for(i = 0; i < nodes; ++i) {
		gx0 = f_state->force[0][i];
		gy0 = f_state->force[1][i];
		rho0 = f_state->rho[i];
		ux0 = f_state->u[0][i];
		uy0 = f_state->u[1][i];

		for(k = 0; k < Q; ++k)
			lbm_state->f[k][i] = lbm_state->f[k][i] -
					omega * (lbm_state->f[k][i] - feq(k, rho0, ux0, uy0)) +
					3.0 * weight[k] * (cx[k]*gx0 + cy[k]*gy0);
	}
}

void stream(FlowState * f_state, LbmState * lbm_state)
{
	unsigned int i, j, lx_m, lx_p, ly_m, ly_p, idx;
	unsigned int nodes = f_state->lx * f_state->ly;

	#pragma omp for private(i, j, lx_m, lx_p, ly_m, ly_p, idx)
	for(idx = 0; idx < nodes; ++idx) {
		i = floor(idx / f_state->ly);
		j = idx - i*f_state->ly;

		/** Permute the indices, corresponding to a periodic boundary condition
			This also prevents segfaults :) */
		lx_m = i==0 				? f_state->lx-1 : i-1;
		lx_p = i==(f_state->lx-1) 	? 0 			: i+1;
		ly_m = j==0 				? f_state->ly-1	: j-1;
		ly_p = j==(f_state->ly-1)	? 0				: j+1;

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

void create_node(Node * node, unsigned int idx, FlowState * f_state, LbmState * lbm_state)
{
	unsigned int i, j, k;

	i = floor(idx / f_state->ly);
	j = idx - i*f_state->ly;

	/* Copy coordinates */
	node->coord[0] = i;
	node->coord[1] = j;

	/* Copy density */
	node->rho = f_state->rho[idx];

	/* Copy force and velocity */
	for(k = 0; k < DIM; ++k) {
		node->u[k] = f_state->u[k][idx];
		node->force[k] = f_state->force[k][idx];
	}

	/* Copy the post-stream distribution functions */
	for(k = 0; k < Q; ++k)
		node->f[k] = lbm_state->f_next[k][idx];
}

void hydrovar_and_bcs(FlowState * f_state, LbmState * lbm_state)
{
	unsigned int k, idx;
	unsigned int nodes = f_state->lx * f_state->ly;

	/* Evaluate the hydrodynamic variables */
	#pragma omp for private(k, idx)
	for(idx = 0; idx < nodes; ++idx) {
		if(f_state->macro_bc[idx] != 0) {
			Node node;
			create_node(&node, idx, f_state, lbm_state);
			macro_bc(&node, f_state, f_state->macro_bc[idx]);

			/* Copy back the state to the global arrays */
			f_state->rho[idx] = node.rho;
			f_state->u[0][idx] = node.u[0];
			f_state->u[1][idx] = node.u[1];
		} else {
			f_state->u[0][idx] = 0.0;
			f_state->u[1][idx] = 0.0;
			f_state->rho[idx] = 0.0;

			for(k = 0; k < Q; ++k) {
				f_state->u[0][idx] += cx[k] * lbm_state->f_next[k][idx];
				f_state->u[1][idx] += cy[k] * lbm_state->f_next[k][idx];
				f_state->rho[idx] += lbm_state->f_next[k][idx];
			}

			f_state->u[0][idx] /= f_state->rho[idx];
			f_state->u[1][idx] /= f_state->rho[idx];
		}

		/* Implement microscopic bcs */
		if(f_state->micro_bc[idx] != 0) {
			Node node;
			create_node(&node, idx, f_state, lbm_state);
			micro_bc(&node, f_state->micro_bc[idx]);

			for(k = 0; k < Q; ++k)
				lbm_state->f_next[k][idx] = node.f[k];
		}
	}
}

void swap_states(LbmState * lbm_state)
{
	unsigned int k;
	double * temp[Q];
	for(k = 0; k < Q; ++k) {
		temp[k] = lbm_state->f[k];
		lbm_state->f[k] = lbm_state->f_next[k];
	}

	for(k = 0; k < Q; ++k)
		lbm_state->f_next[k] = temp[k];
}

void lbm_write_state_binary(FILE * handle, const LbmState * lbm_state)
{
	size_t i;

	write_uint(handle, lbm_state->lx);
	write_uint(handle, lbm_state->ly);

	size_t nodes = lbm_state->lx * lbm_state->ly;

	for(i = 0; i < Q; ++i)
		write_n_doubles(handle, lbm_state->f[i], nodes);

	for(i = 0; i < Q; ++i)
		write_n_doubles(handle, lbm_state->f_next[i], nodes);
}

void lbm_read_state_binary(FILE * handle, LbmState * lbm_state)
{
	size_t i;

	read_uint(handle, &lbm_state->lx);
	read_uint(handle, &lbm_state->ly);

	size_t nodes = lbm_state->lx * lbm_state->ly;

	for(i = 0; i < Q; ++i)
		read_n_doubles(handle, lbm_state->f[i], nodes);
	for(i = 0; i < Q; ++i)
		read_n_doubles(handle, lbm_state->f_next[i], nodes);
}
