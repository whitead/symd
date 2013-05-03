#include "integrate.h"



void integrate_1(double time_step, double* positions, double* velocities,
		 double* forces, double* masses, double* box_size, unsigned int n_dims, 
		 unsigned int n_particles) {

  unsigned int i, j, k;
  for(i = 0; i < n_particles; i++) {
    for(j = 0; j < n_dims; j++) {
      k = i * n_dims + j;
      velocities[k] += time_step * forces[k] / (2 * masses[i]);
      positions[k] += time_step * velocities[k];
      positions[k] = wrap(positions[k], box_size[j]);
    }
  }

}

void integrate_2(double time_step, double* positions, double* velocities,
		 double* forces, double* masses, double* box_size, unsigned int n_dims, 
		 unsigned int n_particles) {

  unsigned int i, j, k;
  for(i = 0; i < n_particles; i++) {
    for(j = 0; j < n_dims; j++) {
      k = i * n_dims + j;
      velocities[k] += time_step * forces[k] / (2 * masses[i]);
    }
  }
}
