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
  const double epsilon;
  const double sigma;
  Nlist_Parameters* nlist;
} Lj_Parameters;
#endif 

/*
 * Calculates forces and rebuilds neighbor list
 */
double gather_forces(void* parameters, double* positions, double* forces, double* masses, 
		     double* box_size, unsigned int n_dims, unsigned int n_particles);

void free_forces(void* parameters);

#ifdef HARMONIC
void* build_harmonic(double k);
#endif

#ifdef LJ
void* build_lj(const double epsilon, const double sigma, Nlist_Parameters* nlist);
#endif 

#endif
