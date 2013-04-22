#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "clbm.h"
#include "micro_bc.h"
#include "macro_bc.h"
#include "plot.h"
#include "fsi.h"
#include "macros.h"
#include "data_types.h"
#include "input.h"
#include <time.h>

void init_flow(FlowParams *, FlowState *);
void destroy_flow(FlowState *);
void post_process(unsigned int, FlowParams *, FlowState *, ParticleState *);
void swap_states(LbmState *);

void print_info(FlowParams * params, FsiParams * fsi_params)
{
	double visc = (1.0/3.0) * (params->tau - 0.5);
	double Re = params->u_max * (params->ly/2) / visc;
	double conf = (fsi_params->a / (params->ly/2));
	double Re_p = Re * conf * conf;
	double St = fsi_params->rho / params->rho * Re_p;
	printf("Re_d = %f\n", Re);
	printf("Re_p = %f\n", Re_p);
	printf("St = %f\n", St);
}

int main(int argc, char ** argv)
{
	unsigned int Nt = 100000, it, output_step = 1000;
	clock_t start, end;
	FlowParams flow_params;
	FsiParams fsi_params;
	FlowState flow_state;
	ParticleState particle_state;
	LbmState lbm_state;

	// Set physical parameters
	parse_input(argc, argv, &flow_params, &fsi_params);

	// Print info
	print_info(&flow_params, &fsi_params);

	// Setup plotting engine
	init_plot();

	// Initialize structs
	fsi_init_state(&fsi_params, &particle_state);
	init_flow(&flow_params, &flow_state);
	lbm_init_state(&flow_state, &lbm_state);

	start = clock();

	for(it = 1; it <= Nt; ++it) {
		if((it % output_step) == 0) {
			end = clock();
			printf("Iteration %d/%d, time/iteration %f, angle=%f, ang_vel=%f\n", it, Nt,
					((float)(end - start))/CLOCKS_PER_SEC/output_step, particle_state.angle, particle_state.ang_vel);
			start = end;
		}

		// Set reference velocity
		flow_state.u_ref = flow_params.u_max * sin(flow_params.f * it);

		// Solve the fsi problem
		fsi_run(&flow_state, &particle_state);

		// Solve the flow problm
		lbm_run(&flow_state, &lbm_state);

		// Post process the result
		post_process(it, &flow_params, &flow_state, &particle_state);

		// Swap the f_next and f array in the LbmState struct
		swap_states(&lbm_state);
	}

	// Clean up
	fsi_destroy_state(&particle_state);
	lbm_destroy_state(&lbm_state);
	destroy_flow(&flow_state);
	destroy_plot();

	return 0;
}

void init_flow(FlowParams * f_params, FlowState * f_state)
{
	unsigned int i, j, k, nx, ny, idx;

	// Macroscopic initial conditions
	nx = f_params->lx;
	ny = f_params->ly;

	f_state->u_ref = f_params->u_max;
	f_state->lx = nx;
	f_state->ly = ny;
	f_state->force[0] = (double *) malloc(nx * ny * sizeof(double));
	f_state->force[1] = (double *) malloc(nx * ny * sizeof(double));
	f_state->u[0] = (double *) malloc(nx * ny * sizeof(double));
	f_state->u[1] = (double *) malloc(nx * ny * sizeof(double));
	f_state->rho = (double *) malloc(nx * ny * sizeof(double));
	f_state->macro_bc = (int *) malloc(nx * ny * sizeof(int));
	f_state->micro_bc = (int *) malloc(nx * ny * sizeof(int));
	f_state->is_corner = (int *) malloc(nx * ny * sizeof(int));
	f_state->tau = f_params->tau;

	// Initialize arrays
	for(i = 0; i < nx; ++i) {
		for(j=0; j < ny; ++j) {
			idx = i*ny + j;

			f_state->u[0][idx] = 0;
			f_state->u[1][idx] = 0;

			for(k = 0; k < DIM; ++k)
				f_state->force[k][idx] = 0;

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

void destroy_flow(FlowState * f_state)
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

void post_process(unsigned int it, FlowParams * f_params, FlowState * f_state, ParticleState * p_state)
{
	if((it % 1000) == 0) {
		PlotOptions plot_opts;
		plot_opts.min_val = -f_params->u_max;
		plot_opts.max_val = f_params->u_max;
		imagesc(f_state->lx, f_state->ly, f_state->u[0], &plot_opts);
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
