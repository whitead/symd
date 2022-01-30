#include "force.h"

typedef struct
{
  double k;
} harmonic_parameters_t;

double h_gather_forces(run_params_t *params, double *positions, double *forces)
{

  unsigned n_dims = params->n_dims;
  unsigned n_particles = params->n_particles;
  force_t *force = params->force_parameters;
  double k = ((harmonic_parameters_t *)force->parameters)->k;
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

void h_free_forces(force_t *force)
{
  harmonic_parameters_t *h_parameters = (harmonic_parameters_t *)force->parameters;
  free(h_parameters);
  free(force);
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
