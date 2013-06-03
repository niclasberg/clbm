#include "lyapunov.h"
#include "fsi.h"
#include <math.h>
#include <stdlib.h>

void lyapunov_init_state(unsigned int it, ParticleState * p_state, LyapunovParticleState * lp_state)
{
	// Initial disturbance and tolerance
	lp_state->d0 = 1.0e-6;
	lp_state->norm_tol = 1.0e-1;

	// Initialize scalar values
	lp_state->cum_sum = 0;
	lp_state->t0 = it;

	// Save a pointer to the base state
	lp_state->base_state = p_state;

	// Copy the state and introduce a small disturbance
	lp_state->perturbed_state = (ParticleState *) malloc(sizeof(ParticleState));
	fsi_copy_state(lp_state->perturbed_state, p_state);
	lp_state->perturbed_state->angle += lp_state->d0;
	fsi_update_particle_nodes(lp_state->perturbed_state);
}

void lyapunov_destroy_state(LyapunovParticleState * lp_state)
{
	fsi_destroy_state(lp_state->perturbed_state);
	free(lp_state->perturbed_state);
	lp_state->base_state = NULL;
}

void lyapunov_run(unsigned int it, LyapunovParticleState * lp_state, FlowState * f_state)
{
	double d, d_square, alpha_square, alpha;

	// Compute torque and update the position of the perturbed state
	fsi_compute_force_on_particle(f_state, lp_state->perturbed_state);
	fsi_update_particle(lp_state->perturbed_state);

	// Find the distance in phase space between the base and the perturbed state
	d = sqrt(pow(lp_state->perturbed_state->angle - lp_state->base_state->angle, 2) +
			   pow((lp_state->perturbed_state->ang_vel/f_state->G - lp_state->base_state->ang_vel/f_state->G), 2));
	alpha = d / lp_state->d0;

	// Normalize if the distance between the trajectories is too great
	if((d < lp_state->d0*lp_state->norm_tol) || d > (lp_state->d0/lp_state->norm_tol)) {
		// Push the perturbed orbit towards the base orbit
		lp_state->perturbed_state->angle = lp_state->base_state->angle + (lp_state->perturbed_state->angle - lp_state->base_state->angle) / alpha;
		lp_state->perturbed_state->ang_vel = lp_state->base_state->ang_vel + (lp_state->perturbed_state->ang_vel*f_state->G - lp_state->base_state->ang_vel*f_state->G) / alpha;

		printf("%e, alpha = %e \n", d, alpha);
		d = sqrt(pow(lp_state->perturbed_state->angle - lp_state->base_state->angle, 2) +
					   pow((lp_state->perturbed_state->ang_vel/f_state->G - lp_state->base_state->ang_vel/f_state->G), 2));
		printf("d after = %e\n", d);

		fsi_update_particle_nodes(lp_state->perturbed_state);

		// Update the lyapunov exponent
		lp_state->cum_sum += log(alpha);
		lp_state->lambda = lp_state->cum_sum / (f_state->G*((double)it - lp_state->t0));

		printf("Lyapunov exponent: %f\n", lp_state->lambda);
	}
}
