#include <gsl/gsl_rng.h>

#ifndef THERMOSTAT_H_
#define THERMOSTAT_H_



/*
 * Returns its conserved quantity. 0 for thermostats which don't coserve anything
 */
double thermostat(double temperature, double time_step, void* thermostat_parameters, double* positions, double* velocities, double* masses, unsigned int n_dims, unsigned int n_particles);

void free_thermostat(void* thermostat_parameters);

#ifdef ANDERSON
#define THERMOSTAT

typedef struct {
  gsl_rng * rng;
  double nu;
} Anderson_Params;

void* build_anderson(unsigned int seed, double collision_freq);

#endif

#ifdef BUSSI
#define THERMOSTAT

typedef struct {
  gsl_rng * rng;
  double taut;
} Bussi_Params;

void* build_bussi(unsigned int seed, double taut);

#endif

#endif //THERMOSTAT_H_
