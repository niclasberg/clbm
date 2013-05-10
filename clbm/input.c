#include "input.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "macros.h"

void parse_input(int argc, char ** argv, FlowParams * flow_params, FsiParams * fsi_params, OutputParams * output_params)
{
	// Default parameters
	InputParameters params;
	params.kb = 0.5;
	params.St = 10;
	params.freq = 0.1;
	params.Re_p = 1;
	params.u_max = 0.01;
	params.lx = 100;
	params.ly = 100;
	params.conf = 0.2;
	params.init_angle = 0;
	params.init_ang_vel = 0;
	params.output_step = 400;
	params.print_particle_state = 1;
	params.print_rho = 0;
	params.print_ux = 0;
	params.print_uy = 0;
	params.timesteps = 12000000;

	// Read from file is the optional argument is supplied
	if(argc > 1 && strlen(argv[1]) > 0)
		read_input_file(argv[1], &params);

	// Compute resulting parameters
	// Flow parameters
	double d = params.ly/2.0; 							// Channel half-height
	double G = (double) params.u_max / d;				// Shear rate (LB units)
	double Re =  params.Re_p / pow(params.conf, 2);		// Channel Reynolds number
	double visc = params.u_max * d / Re;				// Viscocity

	flow_params->tau = 3.0*visc + 0.5;
	flow_params->lx = params.lx;
	flow_params->ly = params.ly;
	flow_params->rho = 1;
	flow_params->f = G * params.freq;
	flow_params->u_max = params.u_max;

	// Fsi parameters
	double a = params.conf * d;
	double b = a * params.kb;

	fsi_params->a = a;
	fsi_params->b = b;
	fsi_params->rho = flow_params->rho * params.St / params.Re_p;
	fsi_params->coord_c[0] = params.lx / 2.0 - 0.5;
	fsi_params->coord_c[1] = params.ly / 2.0 - 0.5;
	fsi_params->init_angle = params.init_angle;
	fsi_params->init_ang_vel = params.init_ang_vel;

	// Compute the number of nodes so that the distance between each
	// node is of order 1
	double h = pow((a - b) / (a + b), 2);
	double circ = PI * (a + b) *
			(1.0 + 3.0*h/(10.0 + sqrt(4.0-3.0*h)));
	fsi_params->nodes = ceil(circ);

	printf("%d", fsi_params->nodes);

	// Output parameters
	output_params->output_step = params.output_step;
	output_params->print_particle_state = params.print_particle_state;
	output_params->print_rho = params.print_rho;
	output_params->print_ux = params.print_ux;
	output_params->print_uy = params.print_uy;
	output_params->timesteps = params.timesteps;
}

void read_input_file(char * file_name, InputParameters * params)
{
	FILE * handle = fopen(file_name, "r");
	if(handle == NULL) {
		fprintf(stderr, "Could not open file %s!\n", file_name);
		exit(1);
	}

	printf("Reading input from %s\n", file_name);

	char line[256];
	char * key, * value;
	char * search = "\t ";

	while(fgets(line, sizeof(line), handle) != NULL) {
		if(strlen(line) < 2)
			continue; // Empty line

		if(line[0] == '#')
			continue;

		// Read key
		key = strtok(line, search);
		value = strtok(NULL, search);

		if(strcmp(key, "Re_p") == 0)
			params->Re_p = atof(value);
		else if(strcmp(key, "conf") == 0)
			params->conf = atof(value);
		else if(strcmp(key, "kb") == 0)
			params->kb = atof(value);
		else if(strcmp(key, "St") == 0)
			params->St = atof(value);
		else if(strcmp(key, "freq") == 0)
			params->freq = atof(value);
		else if(strcmp(key, "umax") == 0)
			params->u_max = atof(value);
		else if(strcmp(key, "freq") == 0)
			params->freq = atof(value);
		else if(strcmp(key, "lx") == 0)
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
		else if(strcmp(key, "timesteps") == 0)
			params->timesteps = atoi(value);
		else {
			fprintf(stderr, "Unknown key %s in file %s\n", key, file_name);
			exit(1);
		}
	}

	fclose(handle);
}
