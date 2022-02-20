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

thermostat_t *build_anderson(double collision_freq, gsl_rng *rng);

thermostat_t *build_bussi(double taut, gsl_rng *rng);

thermostat_t *build_baoab(double gamma, gsl_rng *rng);

#endif // THERMOSTAT_H_
