#include <stdlib.h>
#include <math.h>
#include "min_image.h"
#include "nlist.h"

#ifndef FORCE_H_
#define FORCE_H_

typedef struct
{
  double k;
} harmonic_parameters_t;

typedef struct
{
  const double epsilon;
  const double sigma;
  Nlist_Parameters *nlist;
} Lj_Parameters;

/*
 * Calculates forces and rebuilds neighbor list
 */
typedef double (*gather_forces_t)(void *parameters,
                                  double *positions,
                                  double *forces,
                                  double *masses,
                                  double *box_size,
                                  unsigned int n_dims,
                                  unsigned int n_particles);
//                      double *box_size, unsigned int n_dims, unsigned int n_particles)
// double gather_forces(void *parameters, double *positions, double *forces, double *masses,
//                      double *box_size, unsigned int n_dims, unsigned int n_particles);

typedef void (*free_forces_t)(void *parameters);

typedef struct
{
  gather_forces_t gather;
  free_forces_t free;
  void *parameters;
} force_t;

force_t *build_harmonic(double k);

force_t *build_lj(const double epsilon, const double sigma, Nlist_Parameters *nlist);

#endif
