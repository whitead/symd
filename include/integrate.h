#include "util.h"
#include "min_image.h"

/*
 * The first integration should update positions and velocities with last time step's forces.
 */
void integrate_1(double time_step, double* positions, double* velocities,
		 double* forces, double* masses, double* box_size, unsigned int n_dims, 
		 unsigned int n_particles);

/*
 * Update velocities with this time step's forces.
 */
void integrate_2(double time_step, double* positions, double* velocities,
		 double* forces, double* masses, double* box_size, unsigned int n_dims, 
		 unsigned int n_particles);
