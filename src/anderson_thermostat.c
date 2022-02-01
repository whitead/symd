#include "thermostat.h"
#include <gsl/gsl_randist.h>
#include <math.h>

double a_thermostat(thermostat_t *params, double temperature, double time_step, double *positions, double *velocities, double *masses, unsigned int n_dims, unsigned int n_particles)
{

  unsigned int i, j;
  for (i = 0; i < n_particles; i++)
  {
    //randomly select particle
    if (gsl_ran_flat(params->rng, 0, 1) < params->param * time_step)
    {
      //draw its velocity from maxwell-boltz
      for (j = 0; j < n_dims; j++)
      {
        velocities[i * n_dims + j] = gsl_ran_gaussian_ziggurat(params->rng, sqrt(temperature / masses[i]));
      }
    }
  }

  return 0;
}

void a_free_thermostat(thermostat_t *params)
{
  free(params);
}

thermostat_t *build_anderson(double collision_freq, gsl_rng *rng)
{
  thermostat_t *params = (thermostat_t *)malloc(sizeof(thermostat_t));
  params->param = collision_freq;
  params->rng = rng;
  params->thermo_fxn = &a_thermostat;
  params->free = &a_free_thermostat;
  return params;
}
