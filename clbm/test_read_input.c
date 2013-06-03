#include "input.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char ** argv)
{
	InputParameters * params;
	size_t params_size;
	size_t i;

	read_input_file("parameters.txt", &params, &params_size);

	for(i = 0; i < params_size; ++i) {
		printf("-----\n");
		printf("St = %f\n", params[i].St);
		printf("Re = %f\n", params[i].Re_p);
		printf("f = %f\n", params[i].freq);
		printf("kb = %f\n", params[i].kb);
	}

	free(params);

	return 0;
}
