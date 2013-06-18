#ifndef FLOW_H_
#define FLOW_H_
#include "data_types.h"
#include "macros.h"

void flow_init_state(FlowParams *, FlowState *);
void flow_destroy_state(FlowState *);
void flow_copy_state(FlowState *, FlowState *);

#endif /* FLOW_H_ */
