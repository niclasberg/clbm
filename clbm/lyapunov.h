#ifndef LYAPUNOV_H_
#define LYAPUNOV_H_
#include "data_types.h"

void lyapunov_init_state(unsigned int, ParticleState *, LyapunovParticleState *);
void lyapunov_run(unsigned int, LyapunovParticleState *, FlowState *);
void lyapunov_destroy_state(LyapunovParticleState *);

#endif /* LYAPUNOV_H_ */
