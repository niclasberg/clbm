#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "clbm.h" 		/* for lbm_init_state, lbm_destroy_state, lbm_run */
#include "fsi.h" 		/* for fsi_init_state, fsi_destroy_state, fsi_run */
#include "macros.h"		/* for lattice definitions (Q) */
#include "data_types.h"	/* for solution and parameter structs */
#include "input.h" 		/* for parse_input */
#include "output.h" 	/* for write_output */
#include "workerpool.h"
#include "lyapunov.h"
#include "flow.h"

void flow_init_state(FlowParams *, FlowState *);
void flow_destroy_state(FlowState *);
void swap_states(LbmState *);
void print_info(FlowParams *, FsiParams *);
void solve(void *);

int main(int argc, char ** argv)
{
	size_t i;
	InputParameters * params;
	size_t params_count;

	// read input
	if(argc < 2 || strlen(argv[1]) < 2) {
		fprintf(stderr, "Missing parameter file. Usage: ./lbm parameterfile\n");
		exit(-1);
	}

	read_input_file(argv[1], &params, &params_count);

	// Queue up jobs
	workerpool_init(1);

	for(i = 0; i < params_count; ++i)
		workerpool_push_job(solve, (void *) &params[i]);

	// Run all calculations
	workerpool_run();

	// Clean up
	workerpool_destroy();

	return 0;
}

void solve(void * args) {
	InputParameters * params = (InputParameters *) args;

	int is_done = 0;
	unsigned int it, iterations;
	FlowParams flow_params;
	FsiParams fsi_params;
	OutputParams output_params;
	FlowState flow_state;
	ParticleState particle_state;
	LbmState lbm_state;
	LyapunovParticleState * lya_particle_state = NULL;

	// Parse input parameters
	parse_input(params, &flow_params, &fsi_params, &output_params);

	// Setup output
	init_output(&output_params);

	// Print parameter file
	write_parameters(&output_params, params);

	// Initialize structs
	fsi_init_state(&fsi_params, &particle_state);
	flow_init_state(&flow_params, &flow_state);
	lbm_init_state(&flow_state, &lbm_state);

	// Write initial particle state
	write_output(0, &output_params, &flow_state, &particle_state, lya_particle_state);

	// Iteration at which Lyapunov exponent calculation should start (t = 10 * St)
	//unsigned int lya_start_it = ceil(50.0 * params->alpha * params->Re_p / flow_params.G);
	iterations = ceil((20.0*params->alpha*params->Re_p + 1000*2.0*M_PI/params->freq) / flow_params.G);

	for(it = 1; it <= iterations; ++it) {
		// Set reference velocity (used in the boundary conditions)
		flow_state.u_ref = flow_params.u_max * sin(flow_params.f * it);

		// Solve the fsi problem
		fsi_run(&flow_state, &particle_state);

		// Compute lyapunov exponent
		/*if(it > lya_start_it) {
			if( ! lya_particle_state) {
				lya_particle_state = (LyapunovParticleState *) malloc(sizeof(LyapunovParticleState));
				lyapunov_init_state(it, lya_particle_state, &particle_state, &flow_state, &lbm_state);
			} else {
				lyapunov_run(it, lya_particle_state);
				//printf("Lyapunov exponent: %f\n", lya_particle_state->lambda);
				if(lya_particle_state->norm_count >= 100)
					is_done = 1;
			}
		}*/

		// Solve the flow problem
		lbm_run(&flow_state, &lbm_state);

		// Post process the result
		if((it % output_params.output_step) == 0)
			write_output(it, &output_params, &flow_state, &particle_state, lya_particle_state);

		// Swap the f_next and f array in the LbmState struct
		swap_states(&lbm_state);

		// Update timestep
		//++it;
	}

	// Clean up
	destroy_output(&output_params);
	fsi_destroy_state(&particle_state);
	lbm_destroy_state(&lbm_state);
	flow_destroy_state(&flow_state);

	if(lya_particle_state) {
		lyapunov_destroy_state(lya_particle_state);
		free(lya_particle_state);
	}
}

void print_info(FlowParams * flow_params, FsiParams * fsi_params)
{
	double visc = (1.0/3.0) * (flow_params->tau - 0.5);
	double Re = flow_params->u_max * (flow_params->ly/2) / visc;
	double conf = (fsi_params->a / (flow_params->ly/2));
	double Re_p = Re * conf * conf;
	double St = fsi_params->rho / flow_params->rho * Re_p;
	printf("Domain dimensions = %d x %d\n", flow_params->lx, flow_params->ly);
	printf("Re_d = %f\n", Re);
	printf("Re_p = %f\n", Re_p);
	printf("St = %f\n", St);
	printf("Particle major axis length = %f\n", fsi_params->a);
	printf("Particle minor axis length = %f\n", fsi_params->b);
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
