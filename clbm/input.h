#ifndef INPUT_H_
#define INPUT_H_
#include "data_types.h"
#include <string.h>

void input_parse_input_params(InputParameters *, FlowParams *, FsiParams *, OutputParams *);
void input_read_param_file(char *, InputParameters **, size_t *);
void input_init_params(InputParameters * params);
int input_set_parameter(InputParameters *, char *, char *);

#endif /* INPUT_H_ */
