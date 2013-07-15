#ifndef OUTPUT_H_
#define OUTPUT_H_
#include "data_types.h"

void output_init(OutputParams *);
void output_destroy(OutputParams *);
void output_write_state_to_file(unsigned int, OutputParams *, FlowState *, ParticleState *, LyapunovState *);
void output_write_parameters_to_file(OutputParams *, InputParameters *);
static void write_array_to_new_file(char *, unsigned int, unsigned int, double *);

#endif /* OUTPUT_H_ */
