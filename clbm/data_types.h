#ifndef DATA_TYPES_H_
#define DATA_TYPES_H_
#include "macros.h"

typedef struct {
	double f, tau, u_max, rho;
	unsigned int lx, ly;
} FlowParams;

typedef struct {
	double a, b, rho, init_angle, init_ang_vel;
	double coord_c[DIM];
	unsigned int nodes;
} FsiParams;

typedef struct {
	unsigned int output_step;
	int print_ux;
	int print_uy;
	int print_rho;
	int print_particle_state;
} OutputParams;

typedef struct {
	double * f[Q], * f_next[Q];
} LbmState;

typedef struct {
	unsigned int lx, ly;	// Number of nodes
	int * macro_bc;			// Macroscopic bcs, see bcs.h for more info
	int * micro_bc;			// Microscopic bcs, see bcs.h for more info
	int * is_corner;		// Location of corner nodes
	double * force[DIM];	// Body forces
	double * rho;			// Density
	double * u[DIM];		// Velocity field
	double tau;				// Relaxation parameter
	double u_ref;			// Reference velocity
} FlowState;

typedef struct {
	double angle, ang_vel; 	// Current angle and angular velocity
	double * volume;			// Volume occupied by each node
	double * coord_p[DIM];	// Coordinates in the particle reference frame
	double * coord_a[DIM];	// Coordinates in the fixed frame
	unsigned int nodes; 	// Number of nodes
	double coord_c[DIM]; 		// Center-point coordinate
	double inertia;			// Moment of inertia
	double * force_fsi[DIM];// Fluid-structure interaction force
	double width;			// Width of the particle
} ParticleState;

typedef struct {
	double f[Q];
	double u[DIM];
	double rho;
	int coord[DIM];
	double force[DIM];
} Node;

#endif /* DATA_TYPES_H_ */
