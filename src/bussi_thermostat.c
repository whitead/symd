#include "thermostat.h"
#include <gsl/gsl_randist.h>
#include <math.h>

/*
 * Sampling exactly when less than 6, instead of exactly 1. TODO: Check if this
 * is necessary. DONE, it wasn't so I removed it
 */
double resamplekin_sumnoises(gsl_rng* rng, unsigned int nn){
/*
  returns the sum of n independent gaussian noises squared
   (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
*/
  double rr;
  if(nn == 0) {
    return 0.0;
  } else if(nn % 2 == 0) {
    return 2.0 * gsl_ran_gamma(rng, nn/2, 1.0);
  } else {
    rr= gsl_ran_gaussian_ziggurat(rng, 1);
    return 2.0 * gsl_ran_gamma(rng, (nn-1)/2, 1.0) + rr*rr;
  }
}


/* Equation A7 from Bussi2007
  kk:    present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  sigma: target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  ndeg:  number of degrees of freedom of the atoms to be thermalized
  taut:  relaxation time of the thermostat, in units of 'how often this routine is called'
*/

double resamplekin(gsl_rng* rng, double kk, double sigma, unsigned int ndeg, double taut){
  double factor,rr;
  if(taut>0.1){
    factor=exp(-1.0/taut);
  } else{
    factor=0.0;
  }
  rr = gsl_ran_gaussian_ziggurat(rng, 1.0);

  return kk + (1.0-factor)* (sigma*(resamplekin_sumnoises(rng, ndeg-1)+rr*rr)/ndeg-kk)
              + 2.0*rr*sqrt(kk*sigma/ndeg*(1.0-factor)*factor);
}


double thermostat(double temperature, double time_step, void* thermostat_parameters, double* positions, double* velocities, double* masses, unsigned int n_dims, unsigned int n_particles) {

  Bussi_Params* params = (Bussi_Params*) thermostat_parameters;

  unsigned int i, j, ndeg;
  double kenergy = 0;
  double etemp, scaling_factor, new_kenergy;

  //calculate kinetic energy
  for(i = 0; i < n_particles; i++) {
    etemp = 0;
    for(j = 0; j < n_dims; j++) {
      etemp += velocities[i * n_dims + j] * velocities[i * n_dims + j];
    }
    kenergy += 0.5 * etemp * masses[i];
  }
  
  //Removed COM
  ndeg = (n_particles * n_dims - n_dims);
   
  //Equation A7 from Bussi2007
  new_kenergy = resamplekin(params->rng, kenergy, 0.5 * temperature * ndeg, ndeg, params->taut / time_step);

  //according to gromacs implementation, it can be negative due to rounding(?)
  if(new_kenergy < 0)
    new_kenergy = 0;

#ifdef DEBUG
  printf("Bussi kenergy = %f, desired = %f, sampled = %f (scaling factor = %f)\n", kenergy, 0.5 * temperature * ndeg, new_kenergy, sqrt(new_kenergy / kenergy));
#endif
  
  scaling_factor = sqrt(new_kenergy / kenergy);

  for(i = 0; i < n_particles; i++) {
    for(j = 0; j < n_dims; j++) {
      velocities[i * n_dims + j] *= scaling_factor;
    }
  }

  return new_kenergy - kenergy;

}



void* build_bussi(unsigned int seed, double taut) {
  gsl_rng * rng;
  gsl_rng_env_setup();
  rng = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(rng, seed);
  
  Bussi_Params* params = (Bussi_Params*) malloc(sizeof(Bussi_Params));
  params->taut = taut;
  params->rng = rng;
  return((void*) params);

}


void free_thermostat(void* thermostat_parameters) {
 
  Bussi_Params* params = (Bussi_Params*) thermostat_parameters;
  gsl_rng_free(params->rng);
  free(params);
 
}
