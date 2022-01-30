#include "force.h"

double h_gather_forces(void *parameters, double *positions, double *forces, double *masses,
                       double *box_size, unsigned int n_dims, unsigned int n_particles)
{
  double k = ((harmonic_parameters_t *)parameters)->k;
  unsigned int i, j;
  double penergy = 0;
  double r;

#pragma omp parallel for default(shared) private(r) reduction(+ \
                                                              : penergy)
  for (i = 0; i < n_particles; i++)
  {
    r = 0;
    for (j = 0; j < n_dims; j++)
    {
      r += positions[i * n_dims + j] * positions[i * n_dims + j];
      forces[i * n_dims + j] = -k * positions[i * n_dims + j];
    }
    penergy += r / 2;
  }

  return (penergy);
}

void h_free_forces(void *parameters)
{
  harmonic_parameters_t *h_parameters = (harmonic_parameters_t *)parameters;
  free(h_parameters);
}

force_t *build_harmonic(double k)
{

  harmonic_parameters_t *parameters = (harmonic_parameters_t *)malloc(sizeof(harmonic_parameters_t));
  parameters->k = k;
  force_t *harmonic = (force_t *)malloc(sizeof(force_t));
  harmonic->gather = &h_gather_forces;
  harmonic->free = &h_free_forces;
  harmonic->parameters = parameters;
  return harmonic;
}
