#include <gsl/gsl_rng.h>

#ifndef THERMOSTAT_H_
#define THERMOSTAT_H_

typedef struct thermostat_t thermostat_t;

/*
 * Returns its conserved quantity. 0 for thermostats which don't conserve anything
 */
typedef double (*thermostat_fxn_t)(thermostat_t *thermostat_parameters, double temperature, double time_step, double *positions, double *velocities, double *masses, unsigned int n_dims, unsigned int n_particles);

typedef void (*free_thermostat_t)(thermostat_t *thermostat_parameters);

struct thermostat_t
{
  gsl_rng *rng;
  double param;
  thermostat_fxn_t thermo_fxn;
  free_thermostat_t free;
};

thermostat_t *build_anderson(unsigned int seed, double collision_freq);

thermostat_t *build_bussi(unsigned int seed, double taut);

#endif //THERMOSTAT_H_
