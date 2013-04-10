#include<stdlib.h>

#ifndef FORCE_H_
#define FORCE_H_

typedef struct {
  double k;
} Harmonic_Parameters;

double gather_forces(void* parameters, double* positions, double* forces, double* masses, 
		     unsigned int n_dims, unsigned int n_particles);

void* build_harmonic(double k);

#endif
