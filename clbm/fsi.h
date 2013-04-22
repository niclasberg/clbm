#ifndef FSI_H_
#define FSI_H_
#include "data_types.h"
#include "macros.h"

// Public methods
void fsi_destroy(ParticleState *);
void fsi_init_state(FsiParams *, ParticleState *);
void fsi_run(FlowState *, ParticleState *);
void fsi_print_info();

// Implementation methods
static void generate_particle_volume(ParticleState *);
static void rotate_particle(ParticleState *);
static void generate_particle_initial(FsiParams *, ParticleState *);
static double dirac(double, double);
static void print_particle(ParticleState *);
static double pow2(double);

#endif /* FSI_H_ */
