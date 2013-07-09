#ifndef FLOW_H_
#define FLOW_H_
#include "data_types.h"
#include "macros.h"

FlowState * flow_alloc_state(unsigned int, unsigned int);
void flow_init_state(FlowParams *, FlowState *);
void flow_free_state(FlowState *);
FlowState * flow_clone_state(const FlowState *);
void flow_read_state_unformatted(FILE *, FlowState *);
void flow_write_state_unformatted(FILE *, const FlowState *);

#endif /* FLOW_H_ */
