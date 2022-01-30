#include <stdlib.h>
#include <math.h>
#include "min_image.h"
#include "nlist.h"
#include "util.h"

#ifndef FORCE_H_
#define FORCE_H_

typedef struct force_t force_t;

/*
 * Calculates forces and rebuilds neighbor list
 */
typedef double (*gather_forces_t)(run_params_t *parameters,
                                  double *positions,
                                  double *forces);
//                      double *box_size, unsigned int n_dims, unsigned int n_particles)
// double gather_forces(void *parameters, double *positions, double *forces, double *masses,
//                      double *box_size, unsigned int n_dims, unsigned int n_particles);

typedef void (*free_forces_t)(force_t *parameters);

struct force_t
{
  gather_forces_t gather;
  free_forces_t free;
  void *parameters;
};

force_t *build_harmonic(double k);

force_t *build_lj(const double epsilon, const double sigma, nlist_parameters_t *nlist);

force_t *build_soft();
#endif
