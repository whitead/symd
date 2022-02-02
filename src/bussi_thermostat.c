#include "thermostat.h"
#include <gsl/gsl_randist.h>
#include <math.h>

/*
 * Sampling exactly when less than 6, instead of exactly 1. TODO: Check if this
 * is necessary. DONE, it wasn't so I removed it
 */
static inline double resamplekin_sumnoises(gsl_rng *rng, unsigned int nn)
{
  /*
  returns the sum of n independent gaussian noises squared
   (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
*/
  double rr;
  if (nn == 0)
  {
    return 0.0;
  }
  else if (nn % 2 == 0)
  {
    return 2.0 * gsl_ran_gamma(rng, nn / 2, 1.0);
  }
  else
  {
    rr = gsl_ran_gaussian_ziggurat(rng, 1);
    return 2.0 * gsl_ran_gamma(rng, (nn - 1) / 2, 1.0) + rr * rr;
  }
}

/* Equation A7 from Bussi2007
  kk:    present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  sigma: target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  ndeg:  number of degrees of freedom of the atoms to be thermalized
  taut:  relaxation time of the thermostat, in units of 'how often this routine is called'
*/

static inline double resamplekin(gsl_rng *rng, double kk, double sigma, unsigned int ndeg, double taut)
{
  double factor, rr;
  if (taut > 0.1)
  {
    factor = exp(-1.0 / taut);
  }
  else
  {
    factor = 0.0;
  }
  rr = gsl_ran_gaussian_ziggurat(rng, 1.0);

  return kk + (1.0 - factor) * (sigma * (resamplekin_sumnoises(rng, ndeg - 1) + rr * rr) / ndeg - kk) + 2.0 * rr * sqrt(kk * sigma / ndeg * (1.0 - factor) * factor);
}

double b_thermostat(thermostat_t *params, double temperature, double time_step, double *positions, double *velocities, double *masses, unsigned int n_dims, unsigned int n_particles)
{

  unsigned int i, j, ndeg;
  double kenergy = 0;
  double etemp, scaling_factor, new_kenergy;

  //calculate kinetic energy
#pragma omp parallel for default(shared) private(etemp) reduction(+ \
                                                                  : kenergy)
  for (i = 0; i < n_particles; i++)
  {
    etemp = 0;
    for (j = 0; j < n_dims; j++)
    {
      etemp += velocities[i * n_dims + j] * velocities[i * n_dims + j];
    }
    kenergy += 0.5 * etemp * masses[i];
  }
  //Removed COM
  ndeg = (n_particles * n_dims - n_dims);

  //Equation A7 from Bussi2007
  new_kenergy = resamplekin(params->rng, kenergy, 0.5 * temperature * ndeg, ndeg, params->param / time_step);

  //according to gromacs implementation, it can be negative due to rounding(?)
  if (new_kenergy < 0)
    new_kenergy = 0;

#ifdef DEBUG
  printf("Bussi kenergy = %f, desired = %f, sampled = %f (scaling factor = %f)\n", kenergy, 0.5 * temperature * ndeg, new_kenergy, sqrt(new_kenergy / kenergy));
#endif

  scaling_factor = sqrt(new_kenergy / kenergy);
  if (scaling_factor != scaling_factor)
    scaling_factor = 1;

    //update velocities
#pragma omp parallel for
  for (i = 0; i < n_particles; i++)
  {
    for (j = 0; j < n_dims; j++)
    {
      velocities[i * n_dims + j] *= scaling_factor;
    }
  }

  return new_kenergy - kenergy;
}

void b_free_thermostat(thermostat_t *params)
{
  free(params);
}

thermostat_t *build_bussi(double taut, gsl_rng *rng)
{
  thermostat_t *params = (thermostat_t *)malloc(sizeof(thermostat_t));
  params->param = taut;
  params->rng = rng;
  params->thermo_fxn = &b_thermostat;
  params->free = &b_free_thermostat;
  return params;
}
