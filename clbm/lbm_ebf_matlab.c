#include <mex.h>
#include "flow.h"
#include "fsi.h"
#include "clbm.h"
#include "macros.h"
#include "data_types.h"
#include "input.h"
#include <float.h>
#include <stdlib.h>
#include <math.h>

#define STOP_AFTER_ITERATIONS 0
#define HOMOCLINIC_ORBIT_DETERMINATION 1

void swap_states(LbmState *);
void read_matlab_input(InputParameters *, unsigned int *, unsigned int *, const mxArray *);
void print_params(InputParameters *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	unsigned int it;
	unsigned int pre_it;
	unsigned int mode;

	/* Check input */
	if(nrhs != 1)
		mexErrMsgTxt("Too few input arguments");

	if(! mxIsStruct(prhs[0]))
		mexErrMsgTxt("Expected a struct as input argument");

	/* Translate matlab input to a format that the code understands */
	InputParameters params;
	read_matlab_input(&params, &pre_it, &mode, prhs[0]);
	/* print_params(&params); */

	/* Structs for the parameters */
	FlowParams flow_params;
	FsiParams fsi_params;
	OutputParams output_params;

	/* Parse input parameters */
	input_parse_input_params(&params, &flow_params, &fsi_params, &output_params);

	/* Allocate state structs for the flow and particle */
	FlowState * flow_state = flow_alloc_state(flow_params.lx, flow_params.ly);
	ParticleState * particle_state = fsi_alloc_state(fsi_params.nodes);
	LbmState * lbm_state = lbm_alloc_state(flow_params.lx, flow_params.ly);

	/* Initialize state */
	fsi_init_state(&fsi_params, particle_state);
	flow_init_state(&flow_params, flow_state);
	lbm_init_state(flow_state, lbm_state);

	/* Relax for a few iterations */
	mexPrintf("Generating initial condition, preiterations: %d \n", pre_it);
	for(it = 0; it < pre_it; ++it) {
		/* Run the solver, but do not update particle position nor velocity */
		fsi_run_keep_particle_steady(flow_state, particle_state);
		lbm_run(flow_state, lbm_state);
	}

	/* Export state at first iteration */
	double * output_ang_ptr = NULL, * output_ang_vel_ptr = NULL;
	unsigned int Nt = output_params.timesteps+1;
	if(nlhs > 0) {
		plhs[0] = mxCreateDoubleMatrix(Nt, 1, mxREAL);
		output_ang_ptr = mxGetPr(plhs[0]);
		output_ang_ptr[0] = particle_state->angle;
	}

	if(nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix(Nt, 1, mxREAL);
		output_ang_vel_ptr = mxGetPr(plhs[1]);
		output_ang_vel_ptr[0] = particle_state->ang_vel / flow_state->G;
	}


	/* Main loop */
	mexPrintf("Entering main loop\n");
	int is_done = 0;
	it = 0;
	while( ! is_done) {
		/* Termination condition */
		if(it >= output_params.timesteps)
			break;
		else if(mode == HOMOCLINIC_ORBIT_DETERMINATION) {
			if(particle_state->ang_vel > 0 || fabs(particle_state->angle - fsi_params.init_angle) > PI)
				break;
		}

		it += 1;
		fsi_run(flow_state, particle_state);
		lbm_run(flow_state, lbm_state);

		if(output_ang_ptr)
			output_ang_ptr[it] = particle_state->angle;
		if(output_ang_vel_ptr)
			output_ang_vel_ptr[it] = particle_state->ang_vel / flow_state->G;
	}

	if(output_ang_ptr)
		mxSetM(plhs[0], it+1);
	if(output_ang_vel_ptr)
		mxSetM(plhs[1], it+1);

	fsi_free_state(particle_state);
	flow_free_state(flow_state);
	lbm_free_state(lbm_state);

	mexPrintf("LBM done\n");
}

void read_matlab_input(InputParameters * params, unsigned int * pre_it, unsigned int * mode, const mxArray * minput)
{
	unsigned int i, n_params = mxGetNumberOfFields(minput);
	input_init_params(params);
	*pre_it = 0;
	*mode = STOP_AFTER_ITERATIONS;

	for(i = 0; i < n_params; ++i) {
		const char * key = mxGetFieldNameByNumber(minput, i);
		double value = mxGetScalar(mxGetFieldByNumber(minput, 0, i));

		if(strcmp(key, "Re_p") == 0)
			params->Re_p = value;
		else if(strcmp(key, "kb") == 0)
			params->kb = value;
		else if(strcmp(key, "iterations") == 0)
			params->timesteps = value;
		else if(strcmp(key, "freq") == 0)
			params->freq = value;
		else if(strcmp(key, "alpha") == 0)
			params->alpha = value;
		else if(strcmp(key, "init_angle") == 0)
			params->init_angle = value;
		else if(strcmp(key, "init_ang_vel") == 0)
			params->init_ang_vel = value;
		else if(strcmp(key, "umax") == 0)
			params->u_max = value;
		else if(strcmp(key, "conf") == 0)
			params->conf = value;
		else if(strcmp(key, "lx") == 0)
			params->lx = value;
		else if(strcmp(key, "ly") == 0)
			params->ly = value;
		else if(strcmp(key, "pre_it") == 0)
			*pre_it = value;
		else if(strcmp(key, "mode") == 0) {
			if(value == STOP_AFTER_ITERATIONS || value == HOMOCLINIC_ORBIT_DETERMINATION)
				*mode = value;
			else
				mxErrMsgTxt("Unknown mode");
		} else
			mexErrMsgTxt("Unknown input parameter");
	}
}

void print_params(InputParameters * params)
{
	mexPrintf("Re_p = %f\n", params->Re_p);
	mexPrintf("alpha = %f\n", params->alpha);
	mexPrintf("conf = %f\n", params->conf);
	mexPrintf("freq = %f\n", params->freq);
	mexPrintf("init_ang_vel = %f\n", params->init_ang_vel);
	mexPrintf("init_angle = %f\n", params->init_angle);
	mexPrintf("kb = %f\n", params->kb);
	mexPrintf("lx = %d\n", params->lx);
	mexPrintf("ly = %d\n", params->ly);
	mexPrintf("tau = %f\n", params->tau);
	mexPrintf("iterations = %d\n", params->timesteps);
	mexPrintf("u_max = %f\n", params->u_max);
}
