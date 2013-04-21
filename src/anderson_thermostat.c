#include "thermostat.h"
#include <gsl/gsl_randist.h>
#include <math.h>

void thermostat(double temperature, double time_step, void* thermostat_parameters, double* positions, double* velocities, double* masses, unsigned int n_dims, unsigned int n_particles) {

  Anderson_Params* params = (Anderson_Params*) thermostat_parameters;

  unsigned int i, j;
  for(i = 0; i < n_particles; i++){
    //randomly select particle
    if(gsl_ran_flat(params->rng, 0,1) < params->nu * time_step) {
      //draw its velocity from maxwell-boltz
      for(j = 0; j < n_dims; j++) {
	velocities[i * n_dims + j] = gsl_ran_gaussian_ziggurat(params->rng, sqrt(temperature / masses[i]));
      }
    }
  }

}

void* build_anderson(unsigned int seed, double collision_freq) {
  gsl_rng * rng;
  gsl_rng_env_setup();
  rng = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(rng, seed);
  
  Anderson_Params* params = (Anderson_Params*) malloc(sizeof(Anderson_Params));
  params->nu = collision_freq;
  params->rng = rng;
  return((void*) params);

}
