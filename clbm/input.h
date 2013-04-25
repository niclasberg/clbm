#ifndef INPUT_H_
#define INPUT_H_
#include "data_types.h"

typedef struct {
	double p_length; 		// particle length (along the longest axis)
	double kb;				// particle aspect ratio
	double St;				// Stokes number
	double init_angle;		// Initial orientation of the particle
	double init_ang_vel;	// Initial rotational velocity of particle
	double freq;			// Oscillation frequency of the walls
	double Re_p;			// Particle Reynolds number
	double u_max;			// Maximal wall velocity amplitude
	unsigned int lx, ly;	// Number of grid points for the flow field discretization
	unsigned int nodes;		// Number of nodes for the particle
	unsigned int output_step;
	unsigned int timesteps;
	int print_ux;
	int print_uy;
	int print_rho;
	int print_particle_state;
} InputParameters;

void parse_input(int, char **, FlowParams *, FsiParams *, OutputParams *);

static void read_input_file(char *, InputParameters *);

#endif /* INPUT_H_ */
