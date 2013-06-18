#ifndef CLBM_H_
#define CLBM_H_
#include "data_types.h"

void lbm_init_state(FlowState *, LbmState *);
void lbm_copy_state(LbmState *, LbmState *);
void lbm_run(FlowState *, LbmState *);
void lbm_destroy_state();
void lbm_lattice_info();

static void collide(FlowState * , LbmState *);
static void stream(FlowState * , LbmState *);
static void hydrovar(FlowState * , LbmState *);
static void eval_hydrovar(FlowState *, LbmState *, Node *);
static void implement_bcs(FlowState * , LbmState *);

#endif /* CLBM_H_ */
