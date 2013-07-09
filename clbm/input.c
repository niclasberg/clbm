#include "input.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "macros.h"

void set_parameter(InputParameters *, char *, char *);

void parse_input(InputParameters * params, FlowParams * flow_params, FsiParams * fsi_params, OutputParams * output_params)
{
	// Compute resulting parameters
	// Flow parameters
	double d  = (params->ly - 1.0) / 2.0;				// Channel half-height
	double Re = params->Re_p / pow(params->conf, 2);	// Channel Reynolds number

	if(params->tau == DBL_MAX) {
		// Wall velocity set
		double visc = params->u_max * d / Re;
		flow_params->u_max = params->u_max;
		flow_params->tau = params->tau = 3.0*visc + 0.5;
	} else {
		// Particle relaxation time set
		double visc = (params->tau - 0.5) / 3.0;
		flow_params->u_max = params->u_max = Re * visc / d;
		flow_params->tau = params->tau;
	}

	double G = (double) flow_params->u_max / d;			// Shear rate (LB units)
	flow_params->G = G;
	flow_params->lx = params->lx;
	flow_params->ly = params->ly;
	flow_params->rho = 1;
	flow_params->f = G * params->freq;

	// Fsi parameters
	double a = params->conf * d;
	double b = a * params->kb;

	fsi_params->a = a;
	fsi_params->b = b;
	fsi_params->rho = flow_params->rho * params->alpha;
	fsi_params->coord_c[0] = flow_params->lx / 2.0 - 0.5;
	fsi_params->coord_c[1] = flow_params->ly / 2.0 - 0.5;
	fsi_params->init_angle = params->init_angle;
	fsi_params->init_ang_vel = params->init_ang_vel;

	// Compute the number of nodes so that the distance between each
	// node can be of order 1
	double h = pow((a - b) / (a + b), 2);
	double circ = PI * (a + b) *
			(1.0 + 3.0*h/(10.0 + sqrt(4.0-3.0*h)));
	fsi_params->nodes = ceil(circ);

	// Output parameters
	output_params->output_step = params->output_step;
	output_params->print_particle_state = params->print_particle_state;
	output_params->print_rho = params->print_rho;
	output_params->print_ux = params->print_ux;
	output_params->print_uy = params->print_uy;
	output_params->timesteps = params->timesteps;
	output_params->print_lyapunov = params->print_lyapunov;
	output_params->lyapunov_calc_step = params->lyapunov_calc_step;
	sprintf(output_params->output_folder, "alpha%dconf%.2fNx%dNy%d/kb%.2f/Re%.2f/f%.2f/angle%.2f",
			(int)params->alpha, params->conf, params->lx, params->ly, params->kb, params->Re_p, params->freq, params->init_angle);
}


/*
 * read_input_file
 * Reads an input file which has tab separated key-value-combinations on the form
 * 	key	value	value	value	....
 * 	key2	value	value	value	....
 * on each line.
 * param_array will contain one InputParameters struct for each combination of the
 * parameters provided, e.g. if the file contains
 * 	Re	0.1	0.2
 * 	St	1	2
 * the resulting array will be of length 4, with (Re, St) = {(0.1, 1), (0.1, 2), (0.2, 1), (0.2, 2)}.
 * Default parameter values are provided where needed.
 *
 * Warning: This method is not thread safe!!!
 */
void read_input_file(char * file_name, InputParameters ** param_array, size_t * param_count)
{
	FILE * handle = fopen(file_name, "r");
	if(handle == NULL) {
		fprintf(stderr, "Could not open file %s!\n", file_name);
		exit(-1);
	}
	printf("Reading input from %s\n", file_name);

	char line[256];
	char * key, * value;
	char * search = "\t ";

	InputParameters * params = (InputParameters *) malloc(sizeof(InputParameters));
	size_t params_size = 1;
	size_t params_capacity = 1;
	size_t chunk_size;
	size_t i;

	// Set default parameters for the first parameter set
	params[0].kb = 0.5;
	params[0].alpha = 1;
	params[0].freq = 0.1;
	params[0].Re_p = 1;
	params[0].u_max = 0.01;
	params[0].tau = DBL_MAX;
	params[0].lx = 100;
	params[0].ly = 100;
	params[0].conf = 0.2;
	params[0].init_angle = 0;
	params[0].init_ang_vel = 0;
	params[0].output_step = 400;
	params[0].print_particle_state = 1;
	params[0].print_rho = 0;
	params[0].print_ux = 0;
	params[0].print_uy = 0;
	params[0].timesteps = (unsigned int)-1;
	params[0].print_lyapunov = 0;
	params[0].lyapunov_calc_step = 100;

	// Read file line by line
	while(fgets(line, sizeof(line), handle) != NULL) {
		if(strlen(line) < 2)
			continue; // Empty line

		if(line[0] == '#')
			continue; // comment

		// Read key
		key = strtok(line, search);
		value = strtok(NULL, search);

		// Set parameter for all existing parameter sets
		for(i = 0; i < params_size; ++i) {
			set_parameter(&params[i], key, value);
		}

		chunk_size = params_size;
		value = strtok(NULL, search);

		while(value != NULL) {
			// Check if the array is full. In that case, allocate twice as much memory and copy the old contents.
			if((params_size + chunk_size) > params_capacity) {
				InputParameters * new_params = (InputParameters *) malloc(2 * params_capacity * sizeof(InputParameters));
				for(i = 0; i < params_size; ++i)
					new_params[i] = params[i];
				free(params);

				params = new_params;
				params_capacity = 2*params_capacity;
			}

			// Copy the parameter set to the end of the array
			for(i = 0; i < chunk_size; ++i)
				params[params_size + i] = params[i];

			// Replace the currently considered parameter's value with the one read from the file
			for(i = 0; i < chunk_size; ++i) {
				set_parameter(&params[params_size + i], key, value);
			}

			params_size += chunk_size;

			value = strtok(NULL, search);
		}
	}

	fclose(handle);

	*param_array = params;
	*param_count = params_size;
}

void set_parameter(InputParameters * params, char * key, char * value)
{
	if(strcmp(key, "Re_p") == 0)
		params->Re_p = atof(value);
	else if(strcmp(key, "conf") == 0)
		params->conf = atof(value);
	else if(strcmp(key, "kb") == 0)
		params->kb = atof(value);
	else if(strcmp(key, "alpha") == 0)
		params->alpha = atof(value);
	else if(strcmp(key, "freq") == 0)
		params->freq = atof(value);
	else if(strcmp(key, "umax") == 0) {
		params->u_max = atof(value);
		params->tau = DBL_MAX;
	} else if(strcmp(key, "tau") == 0) {
		params->tau = atof(value);
		params->u_max = DBL_MAX;
	} else if(strcmp(key, "lx") == 0)
		params->lx = atoi(value);
	else if(strcmp(key, "ly") == 0)
		params->ly = atoi(value);
	else if(strcmp(key, "angle") == 0)
		params->init_angle = atof(value);
	else if(strcmp(key, "ang_vel") == 0)
		params->init_ang_vel = atof(value);
	else if(strcmp(key, "output_step") == 0)
		params->output_step = atoi(value);
	else if(strcmp(key, "print_particle") == 0)
		params->print_particle_state = atoi(value);
	else if(strcmp(key, "print_rho") == 0)
		params->print_rho = atoi(value);
	else if(strcmp(key, "print_ux") == 0)
		params->print_ux = atoi(value);
	else if(strcmp(key, "print_uy") == 0)
		params->print_uy = atoi(value);
	else if(strcmp(key, "print_lyapunov") == 0)
		params->print_lyapunov = atoi(value);
	else if(strcmp(key, "lyapunov_calc_step") == 0)
		params->lyapunov_calc_step = atoi(value);
	else if(strcmp(key, "timesteps") == 0)
		params->timesteps = atoi(value);
	else {
		fprintf(stderr, "Input error: Unknown key %s\n", key);
		exit(-1);
	}
}
