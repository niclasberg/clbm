#include "input.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void parse_input(int argc, char ** argv, FlowParams * flow_params, FsiParams * fsi_params)
{
	// Default parameters
	InputParameters params;
	params.p_length = 12;
	params.kb = 2;
	params.St = 10;
	params.freq = 0.1;
	params.Re_p = 1;
	params.u_max = 0.01;
	params.lx = 100;
	params.ly = 100;
	params.nodes = 50;
	params.init_angle = 0;
	params.init_ang_vel = 0;

	// Read from file is the optional argument is supplied
	if(argc > 1 && strlen(argv[1]) > 0)
		read_input_file(argv[1], &params);

	// Compute resulting parameters
	// Flow parameters
	double G = 2.0 * (double) params.u_max / ((double) params.ly);
	double Re = pow((params.ly/2) / params.p_length, 2) * params.Re_p;
	double visc = params.u_max * params.ly/2 / Re;
	flow_params->tau = 3*visc + 0.5;
	flow_params->lx = params.lx;
	flow_params->ly = params.ly;
	flow_params->rho = 1;
	flow_params->f = G * params.freq;
	flow_params->u_max = params.u_max;

	printf("Nodes: %d\n", params.nodes);

	// Fsi parameters
	fsi_params->a = params.p_length;
	fsi_params->b = params.p_length / params.kb;
	fsi_params->rho = flow_params->rho * params.St / params.Re_p;
	fsi_params->nodes = params.nodes;
	fsi_params->coord_c[0] = params.lx / 2.0 - 0.5;
	fsi_params->coord_c[1] = params.ly / 2.0 - 0.5;
	fsi_params->init_angle = params.init_angle;
	fsi_params->init_ang_vel = params.init_ang_vel;
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
		else if(strcmp(key, "p_length") == 0)
			params->p_length = atof(value);
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
		else if(strcmp(key, "nodes") == 0)
			params->nodes = atoi(value);
		else if(strcmp(key, "angle") == 0)
			params->init_angle = atof(value);
		else if(strcmp(key, "ang_vel") == 0)
			params->init_ang_vel = atof(value);
		else {
			fprintf(stderr, "Unknown key %s in file %s\n", key, file_name);
			exit(1);
		}
	}

	fclose(handle);
}
