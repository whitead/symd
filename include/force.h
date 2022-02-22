#include <stdlib.h>
#include <math.h>
#include "box.h"
#include "nlist.h"
#include "params.h"

#ifndef FORCE_H_
#define FORCE_H_

/*
 * Calculates forces and rebuilds neighbor list
 */
typedef double (*gather_forces_t)(run_params_t *parameters,
                                  double *positions,
                                  double *forces);

typedef void (*free_forces_t)(force_t *parameters);

struct force_t
{
  gather_forces_t gather;
  free_forces_t free;
  void *parameters;
};

force_t *build_harmonic(double k);

force_t *build_nlj(const double epsilon, const double sigma, nlist_parameters_t *nlist);

force_t *build_lj(const double epsilon, const double sigma);

force_t *build_soft();

force_t *build_gravity();
#endif
