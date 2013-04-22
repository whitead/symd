#include <gsl/gsl_rng.h>

#ifndef THERMOSTAT_H_
#define THERMOSTAT_H_

#ifdef ANDERSON
typedef struct {
  gsl_rng * rng;
  double nu;
} Anderson_Params;
#endif

#ifdef BUSSI
typedef struct {
  gsl_rng * rng;
  double taut;
} Bussi_Params;
#endif

void thermostat(double temperature, double time_step, void* thermostat_parameters, double* positions, double* velocities, double* masses, unsigned int n_dims, unsigned int n_particles);

#ifdef ANDERSON
#define THERMOSTAT
void* build_anderson(unsigned int seed, double collision_freq);
#endif

#endif
