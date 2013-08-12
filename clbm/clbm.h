#ifndef CLBM_H_
#define CLBM_H_
#include "data_types.h"

LbmState * lbm_alloc_state(unsigned int, unsigned int);
void lbm_free_state(LbmState *);
void lbm_init_state(FlowState *, LbmState *);
LbmState * lbm_clone_state(const LbmState *);
void lbm_copy_state(const LbmState *, LbmState *);
void lbm_run(FlowState *, LbmState *);
void lbm_lattice_info();
void lbm_write_state_binary(FILE *, const LbmState *);
void lbm_read_state_binary(FILE *, LbmState *);


#endif /* CLBM_H_ */
