#include "output.h"
#include <stdio.h>

FILE * particle_output_file;

void init_output(OutputParams * output_params)
{
	particle_output_file = NULL;
	if(output_params->print_particle_state)
		particle_output_file = fopen("particle.txt", "w");
}

void destroy_output()
{
	if(particle_output_file) {
		fclose(particle_output_file);
		particle_output_file = NULL;
	}
}

void write_output(unsigned int it, OutputParams * output_params, FlowState * f_state, ParticleState * p_state)
{
	char file_name[100];
	if(output_params->print_particle_state) {
		fprintf(particle_output_file, "%d\t%.14g\t%.14g\n", it, p_state->angle, p_state->ang_vel);
		fflush(particle_output_file);
	}

	if(output_params->print_rho) {
		sprintf(file_name, "%s%d.txt", "rho", it);
		write_array_to_new_file(file_name, f_state->lx, f_state->ly, f_state->rho);
	}

	if(output_params->print_ux) {
		sprintf(file_name, "%s%d.txt", "ux", it);
		write_array_to_new_file(file_name, f_state->lx, f_state->ly, f_state->u[0]);
	}

	if(output_params->print_uy) {
		sprintf(file_name, "%s%d.txt", "uy", it);
		write_array_to_new_file(file_name, f_state->lx, f_state->ly, f_state->u[1]);
	}
}

void write_array_to_new_file(char * file_name, unsigned int lx, unsigned int ly, double * arr)
{
	unsigned int i, j;
	FILE * handle = fopen(file_name, "w");

	for(i = 0; i < lx; ++i) {
		for(j = 0; j < ly; ++j) {
			fprintf(handle, "%d\t%d\t%.14g\n", i, j, arr[i*ly+j]);
		}
	}
	fclose(handle);
}
