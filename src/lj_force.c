#include "force.h"

double gather_forces(void* parameters, double* positions, double* forces, double* masses, 
		     unsigned int n_dims, unsigned int n_particles) {
  double epsilon = ((Lj_Parameters*) parameters)->epsilon;
  double sigma = ((Lj_Parameters*) parameters)->sigma;
  unsigned int i, j, k;
  double penergy = 0;
  double r, force, diff;
  double force_vector[n_dims];

  //zero forces
  for(i = 0; i < n_particles; i++)
    for(k = 0; k < n_dims; k++)
      forces[i * n_dims + k] = 0;

  //iterate through all pairs
  for(i = 0; i < n_particles; i++) {
    for(j = i+1; j < n_particles; j++) {
      r = 0;

      //distance between particles
      for(k = 0; k < n_dims; k++) {
	diff = positions[j * n_dims + k] - positions[i * n_dims + k];
	r += diff * diff;
	force_vector[k] = diff;
      }
      r = sqrt(r);
      //LJ force and potential
      force = 4 * epsilon * (6 * pow(sigma / r, 7) - 12 * pow(sigma / r, 13));

      for(k = 0; k < n_dims; k++)  {
	forces[i * n_dims + k] += force / r * force_vector[k];
	forces[j * n_dims + k] -= force / r * force_vector[k];
      }

      penergy += 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
    }
  }  
  return(penergy);  
}

void* build_lj(double epsilon, double sigma) {
  Lj_Parameters* parameters = (Lj_Parameters*) malloc(sizeof(Lj_Parameters));
  parameters->epsilon = epsilon;
  parameters->sigma = sigma;
  return (parameters);
}
