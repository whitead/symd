#include "force.h"

double gather_forces(void* parameters, double* positions, double* forces, double* masses, 
		     unsigned int n_dims, unsigned int n_particles) {
  double k = ((Harmonic_Parameters*) parameters)->k;
  unsigned int i, j;
  double penergy = 0;
  double r;

  for(i = 0; i < n_particles; i++) {

    r = 0;
    for(j = 0; j < n_dims; j++) {
      r += positions[i * n_dims + j] * positions[i * n_dims + j];
      forces[i * n_dims + j] = -k * positions[i * n_dims + j];
    }
    penergy += r  / 2;
        
  }

  return(penergy);

}

void* build_harmonic(double k) {
  Harmonic_Parameters* parameters = (Harmonic_Parameters*) malloc(sizeof(Harmonic_Parameters));
  parameters->k = k;
  return (parameters);
}
