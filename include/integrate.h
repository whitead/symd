#include "util.h"
#include "box.h"

void integrate_pos(double time_step, double *positions, double *velocities, unsigned int n_particles);

void integrate_vel(double time_step, double *velocities,
				   double *forces, double *masses,
				   unsigned int n_particles);
