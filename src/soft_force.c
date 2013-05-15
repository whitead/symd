#include "force.h"
#include <math.h>

#define PI 3.14159265


/* A soft cosine potential -- E = cos(pi * r).
 * F = -sin(pi * r) * x / r
 *
 */

double gather_forces(void* parameters, double* positions, double* forces, double* masses, 
		     double* box_size, unsigned int n_dims, unsigned int n_particles) {

  unsigned int i, j, k;
  double penergy = 0;
  double r;
  double force_vector[n_dims];

  //zero forces
  for(i = 0; i < n_particles; i++)
    for(k = 0; k < n_dims; k++)
      forces[i * n_dims + k] = 0;

  for(i = 0; i < n_particles; i++) {
    for(j = i+1; j < n_particles; j++ ){
	r = 0;
	for(k = 0; k < n_dims; k++) {
	  force_vector[k] = min_image_dist(positions[j * n_dims + k] - positions[i * n_dims + k], box_size[k]);
	  r += force_vector[k] * force_vector[k];
	}
	r = sqrt(r);
	if(r < 0.00000001) { //perturb
	  force_vector[0] += 0.01;
	  r = 0.01;
	}
	if(r < 1) {
	  for(k = 0; k < n_dims; k++) {
	    forces[i * n_dims + k] += -sin(PI *  r) * PI * force_vector[k] / r;
	    forces[j * n_dims + k] +=  sin(PI *  r) * PI * force_vector[k] / r;
	  }
	  penergy += cos(PI * r);
	} else
	  penergy += 0;

      }
  }

  return(penergy);
}



void free_forces(void* parameters){
  return;
}
