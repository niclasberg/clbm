#include <mex.h>
#include "flow.h"
#include "fsi.h"
#include "clbm.h"
#include "macros.h"
#include "data_types.h"
#include "input.h"
#include <float.h>

void swap_states(LbmState *);
void read_matlab_input(InputParameters *, const mxArray *);
void print_params(InputParameters *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	unsigned int it;

	/* Check input */
	if(nrhs != 1)
		mexErrMsgTxt("Too few input arguments");

	if(! mxIsStruct(prhs[0]))
		mexErrMsgTxt("Expected a struct as input argument");

	/* Translate matlab input to a format that the code understands */
	InputParameters params;
	read_matlab_input(&params, prhs[0]);
	/* print_params(&params);*/

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

	flow_state->u_ref = flow_params.u_max;

	/* Relax for a few iterations to generate an initial condition */
	for(it = 0; it < 300; ++it) {
		fsi_compute_force_on_particle(flow_state, particle_state);
		fsi_project_force_on_fluid(flow_state, particle_state);
		lbm_run(flow_state, lbm_state);
		swap_states(lbm_state);
	}

	/* Main loop */
	for( ; it < output_params.timesteps; ++it) {
		fsi_run(flow_state, particle_state);
		lbm_run(flow_state, lbm_state);
		swap_states(lbm_state);
	}

	if(nlhs > 0) {
		plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
		double * ptr_to_data = mxGetPr(plhs[0]);
		ptr_to_data[0] = particle_state->angle;
		ptr_to_data[1] = particle_state->ang_vel / flow_state->G;
	}
}

void read_matlab_input(InputParameters * params, const mxArray * minput)
{
	unsigned int i, n_params = mxGetNumberOfFields(minput);
	input_init_params(params);

	for(i = 0; i < n_params; ++i) {
		const char * key = mxGetFieldNameByNumber(minput, i);
		double value = mxGetScalar(mxGetFieldByNumber(minput, 0, i));

		if(strcmp(key, "Re_p") == 0)
			params->Re_p = value;
		else if(strcmp(key, "kb") == 0)
			params->kb = value;
		else if(strcmp(key, "iterations") == 0)
			params->timesteps = value;
		else if(strcmp(key, "alpha") == 0)
			params->alpha = value;
		else if(strcmp(key, "init_angle") == 0)
			params->init_angle = value;
		else if(strcmp(key, "init_ang_vel") == 0)
			params->init_ang_vel = value;
		else if(strcmp(key, "umax") == 0)
			params->u_max = value;
		else if(strcmp(key, "lx") == 0)
			params->lx = value;
		else if(strcmp(key, "ly") == 0)
			params->ly = value;
		else
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
