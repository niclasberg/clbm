#ifndef LYAPUNOV_H_
#define LYAPUNOV_H_
#include "data_types.h"

void lyapunov_init_state(unsigned int, LyapunovState *, ParticleState *, FlowState *, LbmState *);
void lyapunov_run(unsigned int, LyapunovState *, FlowState *);
void lyapunov_destroy_state(LyapunovState *);

#endif /* LYAPUNOV_H_ */
