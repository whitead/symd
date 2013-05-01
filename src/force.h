#include <stdlib.h>
#include <math.h>
#include "min_image.h"

#ifndef FORCE_H_
#define FORCE_H_

#ifdef HARMONIC
typedef struct {
  double k;
} Harmonic_Parameters;
#endif 

#ifdef LJ
#include "nlist.h"
typedef struct {
  double epsilon;
  double sigma;
} Lj_Parameters;
#endif 

double gather_forces(void* parameters, double* positions, double* forces, double* masses, 
		     double* box_size, unsigned int n_dims, unsigned int n_particles);

#ifdef HARMONIC
void* build_harmonic(double k);
#endif

#ifdef LJ
void* build_lj(double epsilon, double sigma);
#endif 

#endif
