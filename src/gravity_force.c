#include "force.h"

double g_gather_forces(run_params_t *params, double *positions, double *forces)
{

  unsigned n_particles = params->n_particles;
  unsigned int i, j, k;
  double penergy = 0, e;
  double r, diff, force_vector[N_DIMS], force;

  for (i = 0; i < n_particles; i++)
  {
    r = 0;
    for (j = i + 1; j < n_particles; j++)
    {
      for (k = 0; k < N_DIMS; k++)
      {
        diff = positions[j * N_DIMS + k] - positions[i * N_DIMS + k];
        r += diff * diff;
        force_vector[k] = diff;
      }

      r = sqrt(r);
      force = 1 / pow(r + 1, 2);
      penergy += e / 2;
      for (k = 0; k < N_DIMS; k++)
      {
        forces[i * N_DIMS + k] += force / r * force_vector[k];
        if (j < params->n_particles)
        {
          forces[j * N_DIMS + k] -= force / r * force_vector[k];
          penergy += e / 2;
        }
      }
    }
  }

  return (penergy);
}

void g_free_forces(force_t *force)
{
  free(force);
}

force_t *build_gravity(double k)
{

  force_t *gravity = (force_t *)malloc(sizeof(force_t));
  gravity->gather = &g_gather_forces;
  gravity->free = &g_free_forces;
  gravity->parameters = NULL;
  return gravity;
}
