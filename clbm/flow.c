#include "flow.h"
#include <stdlib.h>
#include "micro_bc.h" 	/* for microscopic boundary condition definitions */
#include "macro_bc.h"	/* for microscopic boundary condition definitions */

void flow_init_state(FlowParams * f_params, FlowState * f_state)
{
	unsigned int i, j, k, nx, ny, idx;

	// Domain dimensions
	nx = f_params->lx;
	ny = f_params->ly;
	f_state->lx = nx;
	f_state->ly = ny;
	f_state->G = f_params->G;

	// Physical parameters
	f_state->u_ref = f_params->u_max;
	f_state->tau = f_params->tau;

	// Solution arrays
	f_state->force[0] = (double *) malloc(nx * ny * sizeof(double));
	f_state->force[1] = (double *) malloc(nx * ny * sizeof(double));
	f_state->u[0] = (double *) malloc(nx * ny * sizeof(double));
	f_state->u[1] = (double *) malloc(nx * ny * sizeof(double));
	f_state->rho = (double *) malloc(nx * ny * sizeof(double));

	// Boundary condition arrays
	f_state->macro_bc = (int *) malloc(nx * ny * sizeof(int));
	f_state->micro_bc = (int *) malloc(nx * ny * sizeof(int));
	f_state->is_corner = (int *) malloc(nx * ny * sizeof(int));

	// Initialize arrays
	for(i = 0; i < nx; ++i) {
		for(j=0; j < ny; ++j) {
			idx = i*ny + j;

			for(k = 0; k < DIM; ++k) {
				f_state->force[k][idx] = 0;
				f_state->u[k][idx] = 0;
			}

			f_state->rho[idx] = f_params->rho;
			f_state->macro_bc[idx] = 0;
			f_state->micro_bc[idx] = 0;
			f_state->is_corner[idx] = 0;
		}
	}

	// Set boundary conditions
	for(i = 0; i < nx; ++i) {
		f_state->macro_bc[i*ny + ny-1] = bc_north;
		f_state->macro_bc[i*ny + 0   ] = bc_south;

		f_state->micro_bc[i*ny + ny-1] = bc_regularized_north;
		f_state->micro_bc[i*ny + 0   ] = bc_regularized_south;
	}
}

void flow_copy_state(FlowState * dest, FlowState * src)
{
	unsigned int i, j, k, nx, ny, idx;

	// Domain dimensions
	nx = src->lx;
	ny = src->ly;
	dest->lx = nx;
	dest->ly = ny;
	dest->G = src->G;

	// Physical parameters
	dest->u_ref = src->u_ref;
	dest->tau = src->tau;

	// Solution arrays
	dest->force[0] = (double *) malloc(nx * ny * sizeof(double));
	dest->force[1] = (double *) malloc(nx * ny * sizeof(double));
	dest->u[0] = (double *) malloc(nx * ny * sizeof(double));
	dest->u[1] = (double *) malloc(nx * ny * sizeof(double));
	dest->rho = (double *) malloc(nx * ny * sizeof(double));

	// Boundary condition arrays
	dest->macro_bc = (int *) malloc(nx * ny * sizeof(int));
	dest->micro_bc = (int *) malloc(nx * ny * sizeof(int));
	dest->is_corner = (int *) malloc(nx * ny * sizeof(int));

	// Copy data
	for(i = 0; i < nx*ny; ++i) {
		dest->macro_bc[i] = src->macro_bc[i];
		dest->micro_bc[i] = src->micro_bc[i];
		dest->is_corner[i] = src->is_corner[i];
		dest->rho[i] = src->rho[i];

		for(k = 0; k < DIM; ++k) {
			dest->force[k][i] = src->force[k][i];
			dest->u[k][i] = src->u[k][i];
		}
	}
}

void flow_destroy_state(FlowState * f_state)
{
	unsigned int k;
	for(k = 0; k < DIM; ++k) {
		free(f_state->force[k]);
		free(f_state->u[k]);
	}
	free(f_state->rho);
	free(f_state->macro_bc);
	free(f_state->micro_bc);
	free(f_state->is_corner);
}
