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
#include "flow.h"

void print_info(FlowParams *, FsiParams *);
void solve(void *);
int lbm_ebf_step(unsigned int, FlowParams *, FlowState *, ParticleState *, LbmState *);

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

	input_read_param_file(argv[1], &params, &params_count);

	// Initialize the worker pool
	workerpool_init(1);

	// Queue up jobs
	for(i = 0; i < params_count; ++i)
		workerpool_push_job(solve, (void *) &params[i]);

	// Run all calculations
	workerpool_run();

	// Clean up
	workerpool_destroy();
	free(params);

	return 0;
}

void solve(void * args) {
	unsigned int it, iterations, lya_start_it;
	char checkpoint_file_name[512];

	// Structs for the parameters
	InputParameters * params = (InputParameters *) args;
	FlowParams flow_params;
	FsiParams fsi_params;
	OutputParams output_params;

	// Parse input parameters
	input_parse_input_params(params, &flow_params, &fsi_params, &output_params);

	// Allocate state structs for the flow and particle
	FlowState * flow_state = flow_alloc_state(flow_params.lx, flow_params.ly);
	ParticleState * particle_state = fsi_alloc_state(fsi_params.nodes);
	LbmState * lbm_state = lbm_alloc_state(flow_params.lx, flow_params.ly);

	// State structs for Lyapunov calculation
	LyapunovState * lya_state = NULL;
	FlowState * lya_flow_state = NULL;
	ParticleState * lya_particle_state = NULL;
	LbmState * lya_lbm_state = NULL;

	// Setup output
	output_init(&output_params);

	// Try to initialize the solution from a checkpoint file
	sprintf(checkpoint_file_name, "%s/checkpoint.dat", output_params.output_folder);
	FILE * cp_handle = fopen(checkpoint_file_name, "r");

	if(cp_handle) {
		// Checkpoint file existed, read the file and initialize the states
		int lya_exists;

		read_uint(cp_handle, &it);
		fsi_read_state_binary(cp_handle, particle_state);
		flow_read_state_unformatted(cp_handle, flow_state);
		lbm_read_state_binary(cp_handle, lbm_state);

		read_uint(cp_handle, &lya_exists);
		if(lya_exists) {
			lya_particle_state = fsi_alloc_state(fsi_params.nodes);
			fsi_read_state_binary(cp_handle, lya_particle_state);
			lya_flow_state = flow_alloc_state(flow_params.lx, flow_params.ly);
			flow_read_state_unformatted(cp_handle, lya_flow_state);
			lya_lbm_state = lbm_alloc_state(flow_params.lx, flow_params.ly);
			lbm_read_state_binary(cp_handle, lya_lbm_state);

			lya_state = malloc(sizeof(LyapunovState));
			read_double(cp_handle, &lya_state->d0);
			read_double(cp_handle, &lya_state->cum_sum);
			read_double(cp_handle, &lya_state->lambda);
			read_uint(cp_handle, &lya_state->t0);
		}

		fclose(cp_handle);
	} else {
		// No checkpoint file, initialize normally

		// Print parameter file
		output_write_parameters_to_file(&output_params, params);

		// Initialize from base state
		it = 0;
		fsi_init_state(&fsi_params, particle_state);
		flow_init_state(&flow_params, flow_state);
		lbm_init_state(flow_state, lbm_state);

		// Print initial state
		output_write_state_to_file(0, &output_params, flow_state, particle_state, lya_state);
	}

	// Number of iterations
	lya_start_it = ceil(20.0*params->alpha*params->Re_p / flow_params.G);
	iterations = lya_start_it + ceil(40.0*params->alpha*params->Re_p / flow_params.G);
	if(iterations > (it + output_params.timesteps))
		iterations = (it + output_params.timesteps);

	// Main loop
	while(it < iterations) {
		// Advance the solution
		it += 1;
		if( ! lbm_ebf_step(it, &flow_params, flow_state, particle_state, lbm_state))
			break;

		// Check if the Lyapunov exponent should be calculated
		if(output_params.print_lyapunov && (it >= lya_start_it)) {
			double d, alpha;

			// Initialize Lyapunov calculation stuff, if not already done
			if( ! lya_state) {
				lya_state = malloc(sizeof(LyapunovState));
				lya_state->d0 = 1.0e-6;
				lya_state->cum_sum = 0;
				lya_state->t0 = it;
				lya_state->lambda = 0;

				// Copy states
				lya_lbm_state = lbm_clone_state(lbm_state);
				lya_flow_state = flow_clone_state(flow_state);
				lya_particle_state = fsi_clone_state(particle_state);

				// Perturb the particle state
				lya_particle_state->angle += lya_state->d0;
				fsi_update_particle_nodes(lya_particle_state);
			} else {
				// Advance the perturbed state
				if( ! lbm_ebf_step(it, &flow_params, lya_flow_state, lya_particle_state, lya_lbm_state))
					break;

				// Normalize and compute Lyapunov exponent
				if(((it - lya_state->t0) % output_params.lyapunov_calc_step) == 0) {
					// Calculate the distance between the original and perturbed orbit
					double ang_vel = particle_state->ang_vel / flow_params.G;
					double lya_ang_vel = lya_particle_state->ang_vel / flow_params.G;

					d = sqrt(pow(lya_particle_state->angle - particle_state->angle, 2) + pow((lya_ang_vel - ang_vel), 2));
					alpha = d / lya_state->d0;

					// Push the perturbed orbit towards the base orbit
					lya_particle_state->angle = particle_state->angle + (lya_particle_state->angle - particle_state->angle) / alpha;
					lya_particle_state->ang_vel = flow_params.G * (ang_vel + (lya_ang_vel - ang_vel) / alpha);

					fsi_update_particle_nodes(lya_particle_state);

					// Update the lyapunov exponent
					lya_state->cum_sum += log(alpha);
					lya_state->lambda = lya_state->cum_sum / (lya_flow_state->G*(it - lya_state->t0));
				}
			}
		}

		// Post process the result
		if((it % output_params.output_step) == 0)
			output_write_state_to_file(it, &output_params, flow_state, particle_state, lya_state);
	}

	// Write a checkpoint file so the simulation can be resumed at later times
	cp_handle = fopen(checkpoint_file_name, "w");
	write_uint(cp_handle, it);
	fsi_write_state_binary(cp_handle, particle_state);
	flow_write_state_unformatted(cp_handle, flow_state);
	lbm_write_state_binary(cp_handle, lbm_state);
	if( ! lya_state) {
		write_uint(cp_handle, 0);
	} else {
		write_uint(cp_handle, 1);
		fsi_write_state_binary(cp_handle, lya_particle_state);
		flow_write_state_unformatted(cp_handle, lya_flow_state);
		lbm_write_state_binary(cp_handle, lya_lbm_state);

		write_double(cp_handle, lya_state->d0);
		write_double(cp_handle, lya_state->cum_sum);
		write_double(cp_handle, lya_state->lambda);
		write_uint(cp_handle, lya_state->t0);
	}
	fclose(cp_handle);

	// Clean up
	output_destroy(&output_params);
	fsi_free_state(particle_state);
	lbm_free_state(lbm_state);
	flow_free_state(flow_state);

	if(lya_state) {
		fsi_free_state(lya_particle_state);
		lbm_free_state(lya_lbm_state);
		flow_free_state(lya_flow_state);
		free(lya_state);
	}
}

int lbm_ebf_step(unsigned int it, FlowParams * flow_params, FlowState * flow_state, ParticleState * particle_state, LbmState * lbm_state)
{
	// Set reference velocity (used in the boundary conditions)
	if(flow_params->f < 1.0e-8)
		flow_state->u_ref = flow_params->u_max;
	else
		flow_state->u_ref = flow_params->u_max * sin(flow_params->f * it);

	// Solve the fsi problem
	fsi_run(flow_state, particle_state);

	// Solve the flow problem
	lbm_run(flow_state, lbm_state);

	// Check for NaN (this check relies on a feature of the floating point system, i.e. NaN != NaN always)
	if(particle_state->ang_vel != particle_state->ang_vel) {
		fprintf(stderr, "ERROR: Solution has diverged.\n");
		return 0;
	}

	return 1;
}

void print_info(FlowParams * flow_params, FsiParams * fsi_params)
{
	double visc = (1.0/3.0) * (flow_params->tau - 0.5);
	double Re = flow_params->u_max * ((flow_params->ly - 1.0)/2.0) / visc;
	double conf = (fsi_params->a / ((flow_params->ly - 1.0)/2.0));
	double Re_p = Re * conf * conf;
	double St = fsi_params->rho / flow_params->rho * Re_p;
	printf("Domain dimensions = %d x %d\n", flow_params->lx, flow_params->ly);
	printf("Confinement = %f\n", conf);
	printf("tau = %f\n", flow_params->tau);
	printf("Re_d = %f\n", Re);
	printf("Re_p = %f\n", Re_p);
	printf("St = %f\n", St);
	printf("Particle major axis length = %f\n", fsi_params->a);
	printf("Particle minor axis length = %f\n", fsi_params->b);
}
