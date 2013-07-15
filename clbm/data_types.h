#ifndef DATA_TYPES_H_
#define DATA_TYPES_H_
#include "macros.h"
#include <stdio.h>

typedef struct {
	double conf; 			/* Confinement: particle_major_length / domain_height */
	double kb;				/* Particle aspect ratio */
	double alpha;			/* Density ratio (particle / fluid) */
	double init_angle;		/* Initial orientation of the particle */
	double init_ang_vel;	/* Initial rotational velocity of particle */
	double freq;			/* Oscillation frequency of the walls */
	double Re_p;			/* Particle Reynolds number */
	double u_max;			/* Maximal wall velocity amplitude */
	double tau;				/* LBM relaxation time */
	unsigned int lx, ly;	/* Number of grid points for the flow field discretization */
	unsigned int output_step;
	unsigned int timesteps;
	int print_ux;
	int print_uy;
	int print_rho;
	int print_particle_state;
	int print_lyapunov;
	unsigned int lyapunov_calc_step;
} InputParameters;

typedef struct {
	double G, f, tau, u_max, rho;
	unsigned int lx, ly;
} FlowParams;

typedef struct {
	double a, b, rho, init_angle, init_ang_vel;
	double coord_c[DIM];
	unsigned int nodes;
} FsiParams;

typedef struct {
	unsigned int timesteps;
	unsigned int output_step;
	int print_ux;
	int print_uy;
	int print_rho;
	int print_particle_state;
	int print_lyapunov;
	unsigned int lyapunov_calc_step;
	char output_folder[256];
	FILE * output_file;
} OutputParams;

typedef struct {
	unsigned int lx, ly;
	double * f[Q], * f_next[Q];
} LbmState;

typedef struct {
	unsigned int lx, ly;	/* Number of nodes */
	int * macro_bc;			/* Macroscopic bcs, see bcs.h for more info */
	int * micro_bc;			/* Microscopic bcs, see bcs.h for more info */
	int * is_corner;		/* Location of corner nodes */
	double * force[DIM];	/* Body forces */
	double * rho;			/* Density */
	double * u[DIM];		/* Velocity field */
	double tau;				/* Relaxation parameter */
	double u_ref;			/* Reference velocity */
	double G;				/* Time scale */
} FlowState;

typedef struct {
	double angle, ang_vel; 	/* Current angle and angular velocity */
	double * volume;		/* Volume occupied by each node */
	double * coord_p[DIM];	/* Coordinates in the particle reference frame */
	double * coord_a[DIM];	/* Coordinates in the fixed frame */
	unsigned int nodes; 	/* Number of nodes */
	double coord_c[DIM]; 	/* Center-point coordinate */
	double inertia;			/* Moment of inertia */
	double torque;			/* Torque on particle */
	double * force_fsi[DIM];/* Fluid-structure interaction force */
	double width;			/* Width of the particle */
} ParticleState;

typedef struct {
	unsigned int t0; 		/* Iteration where the lyapunov calculation started */
	double d0; 				/* Magnitude of the initial disturbance */
	double lambda; 			/* Lyapunov exponent */
	double cum_sum; 		/* = "lyapunov exponent" * (t - t0) */
} LyapunovState;

typedef struct {
	double f[Q];
	double u[DIM];
	double rho;
	int coord[DIM];
	double force[DIM];
} Node;

#endif /* DATA_TYPES_H_ */
