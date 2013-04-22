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
#include <time.h>

#define NX 100
#define NY 100

unsigned int index_2d(unsigned int i, unsigned int j)
{
	return i*NY + j;
}

void parse_input(int argc, char **, FlowParams *, FsiParams *);
void init_flow(FlowParams *, FlowState *);
void post_process(unsigned int, FlowState *, ParticleState *);
void swap_states(LbmState *);

void print_info(FlowParams * params, FsiParams * fsi_params)
{
	double visc = (1.0/3.0) * (params->tau - 0.5);
	double Re = params->u_max * NY/2 / visc;
	double conf = (fsi_params->a / (NY/2));
	double Re_p = Re * conf * conf;
	double St = fsi_params->rho / params->rho * Re_p;
	printf("Channel Re = %f\n", Re);
	printf("Particle Re = %f\n", Re_p);
	printf("St = %f\n", St);
}

int main(int argc, char ** argv)
{
	unsigned int Nt = 10000, it, output_step = 1000;
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
	init_flow(&flow_params, &flow_state);
	fsi_init(&fsi_params, &particle_state);
	lbm_init(&flow_state, &lbm_state);
	//bc_init(INIT_DIFFUSIVE, &lbm_state, NX, NY);

	double G = 2.0 * (double) flow_params.u_max / ((double) NY);
	start = clock();

	for(it = 1; it <= Nt; ++it) {
		if((it % output_step) == 0) {
			end = clock();
			printf("Iteration %d/%d, time/iteration %f, angle=%f, ang_vel=%f\n", it, Nt,
					((float)(end - start))/CLOCKS_PER_SEC/output_step, particle_state.angle, particle_state.ang_vel);
			start = end;
		}

		// Set reference velocity
		flow_state.u_ref = flow_params.u_max * sin(G * flow_params.f * it);

		// Solve the fsi problem
		fsi_run(&flow_state, &particle_state);

		// Solve the flow problm
		lbm_run(&flow_state, &lbm_state);

		// Post process the result
		//post_process(it, &flow_state, &particle_state);

		// Swap the f_next and f array in the LbmState struct
		swap_states(&lbm_state);
	}

	// Clean up
	fsi_destroy(&particle_state);
	lbm_destroy(&lbm_state);
	destroy_plot();

	return 0;
}

void parse_input(int argc, char ** argv, FlowParams * flow_params, FsiParams * fsi_params)
{
	// Flow parameters
	flow_params->f = 0.1;
	flow_params->tau = 0.6;
	flow_params->u_max = 0.01;
	flow_params->rho = 1;

	// Fsi parameters
	double a = 12;
	double b = 6;
	fsi_params->a = a;
	fsi_params->b = b;
	fsi_params->rho = 1;
	fsi_params->coord_c[0] = NX / 2.0 - 0.5;
	fsi_params->coord_c[1] = NY / 2.0 - 0.5;
	fsi_params->nodes = 60;

	double G = 2.0 * flow_params->u_max / ((double) NY);
	fsi_params->init_angle = PI / 2.0;
	fsi_params->init_ang_vel = -(G / (a*a + b*b)) *
			((b*cos(fsi_params->init_angle))*(b*cos(fsi_params->init_angle)) +
			 (a*sin(fsi_params->init_angle))*(a*sin(fsi_params->init_angle)));
}

void init_flow(FlowParams * f_params, FlowState * f_state)
{
	unsigned int i, j, k, idx;

	// Macroscopic initial conditions
	f_state->u_ref = f_params->u_max;
	f_state->lx = NX;
	f_state->ly = NY;
	f_state->force[0] = (double *) malloc(NX * NY * sizeof(double));
	f_state->force[1] = (double *) malloc(NX * NY * sizeof(double));
	f_state->u[0] = (double *) malloc(NX * NY * sizeof(double));
	f_state->u[1] = (double *) malloc(NX * NY * sizeof(double));
	f_state->rho = (double *) malloc(NX * NY * sizeof(double));
	f_state->macro_bc = (int *) malloc(NX * NY * sizeof(int));
	f_state->micro_bc = (int *) malloc(NX * NY * sizeof(int));
	f_state->is_corner = (int *) malloc(NX * NY * sizeof(int));
	f_state->tau = f_params->tau;

	// Initialize arrays
	for(i = 0; i < NX; ++i) {
		for(j=0; j < NY; ++j) {
			idx = index_2d(i, j);

			f_state->u[0][idx] = 0;//f_state->u_ref * (-1 + 2 *(double)j / (f_state->ly - 1.0));
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
	for(i = 0; i < (NX); ++i) {
		f_state->macro_bc[index_2d(i, NY-1)] = bc_north;
		f_state->macro_bc[index_2d(i, 0)] = bc_south;

		f_state->micro_bc[index_2d(i, NY-1)] = bc_regularized_north;
		f_state->micro_bc[index_2d(i, 0)] = bc_regularized_south;
	}

	/*for(i = 1; i < (NY-1); ++i) {
		f_state->macro_bc[index_2d(0, i)] = bc_west;
		f_state->macro_bc[index_2d(NX-1, i)] = bc_east;

		f_state->micro_bc[index_2d(0, i)] = bc_regularized_west;
		f_state->micro_bc[index_2d(NX-1, i)] = bc_regularized_east;

	}

	// Corner bcs
	f_state->macro_bc[index_2d(0, 0)] = bc_south_west;
	f_state->macro_bc[index_2d(NX-1, 0)] = bc_south_east;
	f_state->macro_bc[index_2d(NX-1, NY-1)] = bc_north_east;
	f_state->macro_bc[index_2d(0, NY-1)] = bc_north_west;

	f_state->is_corner[index_2d(0, 0)] = 1;
	f_state->is_corner[index_2d(NX-1, 0)] = 1;
	f_state->is_corner[index_2d(NX-1, NY-1)] = 1;
	f_state->is_corner[index_2d(0, NY-1)] = 1;

	f_state->micro_bc[index_2d(0, 0)] = bc_zou_he_south_west;
	f_state->micro_bc[index_2d(NX-1, 0)] = bc_zou_he_south_east;
	f_state->micro_bc[index_2d(NX-1, NY-1)] = bc_zou_he_north_east;
	f_state->micro_bc[index_2d(0, NY-1)] = bc_zou_he_north_west;*/

	/*for(j = NY; j > 0; ) {
		--j;
		for(i = 0; i < NX; ++i) {
			printf("%d, ", f_state->macro_bc[i*NY +j]);
		}
		printf("\n");
	}*/
}

void post_process(unsigned int it, FlowState * f_state, ParticleState * p_state)
{
	if((it % 1000) == 0) {
		PlotOptions plot_opts;
		plot_opts.min_val = -0.01;
		plot_opts.max_val = 0.01;
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
