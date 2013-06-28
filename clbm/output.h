#ifndef OUTPUT_H_
#define OUTPUT_H_
#include "data_types.h"

void init_output(OutputParams *);
void destroy_output(OutputParams *);
void write_output(unsigned int, OutputParams *, FlowState *, ParticleState *, LyapunovState *);
static void write_array_to_new_file(char *, unsigned int, unsigned int, double *);

#endif /* OUTPUT_H_ */
