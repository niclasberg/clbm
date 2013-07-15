#ifndef FSI_H_
#define FSI_H_
#include "data_types.h"
#include "macros.h"

/* Public methods */
ParticleState * fsi_alloc_state(unsigned int);
void fsi_free_state(ParticleState *);
void fsi_init_state(FsiParams *, ParticleState *);
void fsi_run(FlowState *, ParticleState *);
void fsi_print_info(ParticleState *);
ParticleState * fsi_clone_state(const ParticleState *);
void fsi_update_particle_nodes(ParticleState *);
void fsi_compute_force_on_particle(FlowState *, ParticleState *);
void fsi_update_particle(ParticleState *);
void fsi_project_force_on_fluid(FlowState *, ParticleState *);
void fsi_write_state_binary(FILE *, const ParticleState *);
void fsi_read_state_binary(FILE *, ParticleState *);

/* Implementation methods */
static void generate_particle_volume(ParticleState *);
static void generate_particle_initial(FsiParams *, ParticleState *);
static double dirac(double, double);
static void print_particle(ParticleState *);
static double pow2(double);

#endif /* FSI_H_ */
