#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "clbm.h" 		/* for lbm_init_state, lbm_destroy_state, lbm_run */
#include "fsi.h" 		/* for fsi_init_state, fsi_destroy_state, fsi_run */
#include "macros.h"		/* for lattice definitions (Q) */
#include "data_types.h"	/* for solution and parameter structs */
#include "input.h" 		/* for parse_input */
#include "output.h" 	/* for write_output */
#include "workerpool.h"
#include "flow.h"
#include <omp.h>

#define IT_COUNT 50000

int main(int argc, char * argv[])
{
	unsigned int it, j, k;
	InputParameters inp;

	input_init_params(&inp);

	FlowParams flow_params;
	FsiParams fsi_params;
	OutputParams output_params;

	unsigned int num_threads_array[9] = {1, 2, 4, 6, 8, 10, 12, 14, 16};

	for(k = 0; k < 9; ++k) {
		unsigned int num_threads = num_threads_array[k];
		omp_set_num_threads(num_threads);

		printf("#threads = %d\n", num_threads);
		inp.lx = 64;
		inp.ly = 64;

		for(j = 0; j < 3; ++j) {
			input_parse_input_params(&inp, &flow_params, &fsi_params, &output_params);

			FlowState * flow_state = flow_alloc_state(flow_params.lx, flow_params.ly);
			ParticleState * particle_state = fsi_alloc_state(fsi_params.nodes);
			LbmState * lbm_state = lbm_alloc_state(flow_params.lx, flow_params.ly);

			fsi_init_state(&fsi_params, particle_state);
			flow_init_state(&flow_params, flow_state);
			lbm_init_state(flow_state, lbm_state);

			time_t start_time = time(NULL);

			for(it = 0; it < IT_COUNT; ++it) {
				fsi_run(flow_state, particle_state);
				lbm_run(flow_state, lbm_state);
			}

			time_t end_time = time(NULL);

			double time_elapsed = difftime(end_time, start_time);

			printf("Nx = %d, Ny = %d, total time = %f, time/iteration = %f\n", lbm_state->lx, lbm_state->ly, time_elapsed, time_elapsed / IT_COUNT);

			flow_free_state(flow_state);
			fsi_free_state(particle_state);
			lbm_free_state(lbm_state);

			inp.lx *= 2;
			inp.ly *= 2;
		}
	}
}
