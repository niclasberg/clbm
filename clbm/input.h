#ifndef INPUT_H_
#define INPUT_H_
#include "data_types.h"
#include <string.h>

void parse_input(InputParameters *, FlowParams *, FsiParams *, OutputParams *);
void read_input_file(char *, InputParameters **, size_t *);

#endif /* INPUT_H_ */
