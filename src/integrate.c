#include "integrate.h"

void integrate_pos(double time_step, double *positions, double *velocities, unsigned int n_particles)
{

  unsigned int i, j, k;
#pragma omp parallel for default(shared) private(i, j, k)
  for (i = 0; i < n_particles; i++)
  {
    for (j = 0; j < N_DIMS; j++)
    {
      k = i * N_DIMS + j;
      positions[k] += time_step * velocities[k] / 2;
    }
  }
}

void integrate_vel(double time_step, double *velocities,
                   double *forces, double *masses,
                   unsigned int n_particles)
{

  unsigned int i, j, k;
#pragma omp parallel for default(shared) private(i, j, k)
  for (i = 0; i < n_particles; i++)
  {
    for (j = 0; j < N_DIMS; j++)
    {
      k = i * N_DIMS + j;
      velocities[k] += time_step * forces[k] / (2 * masses[i]);
    }
  }
}
