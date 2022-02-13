#include "force.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

typedef struct
{
  const double epsilon;
  const double sigma;
  nlist_parameters_t *nlist;
} lj_parameters_t;

static inline double lj(double r, double epsilon, double sigma, double *energy)
{
  double ri6 = pow(sigma / r, 6);
  double ri12 = ri6 * ri6;
  *energy += 4 * epsilon * (ri6 - ri12);
  return 4 * epsilon * (6 * ri6 * sigma / r - 12 * ri12 * sigma / r);
}

double lj_gather_forces(run_params_t *params, double *positions, double *forces)
{
  unsigned int n_dims = N_DIMS;
  unsigned int n_particles = params->n_particles;
  double *box_size = params->box->box_size;
  force_t *force_p = params->force_parameters;
  lj_parameters_t *parameters = (lj_parameters_t *)force_p->parameters;
  const double epsilon = parameters->epsilon;
  const double sigma = parameters->sigma;
  nlist_parameters_t *nlist = parameters->nlist;

  // update neighbor list
  // TODO: Turn back on
  // update_nlist(positions, box_size, n_dims, n_particles, params->n_ghost_particles, nlist);

  unsigned int i, j, k, n;
  int offset;
  double penergy = 0;
  double r, force, diff;
  double force_vector[n_dims];
  double rcut = sqrt(nlist->rcut);
  double e_shift, e;
  double lj_shift = lj(rcut, epsilon, sigma, &e_shift);

#ifdef DEBUG
  // check_nlist(params, nlist, positions, rcut);
#endif

// zero forces
#pragma omp parallel for
  for (i = 0; i < n_particles; i++)
    for (k = 0; k < n_dims; k++)
      forces[i * n_dims + k] = 0;

      // iterate through all particles

      // This seems strange at first,
      // but it really just distributes
      // the work across the threads
      // The only trick is that the neighborlists
      // are not easy to spread out across the workers,
      // hence the conditionals.

#pragma omp parallel default(shared) private(offset, n, i, j, k, r, force, e, force_vector, diff) \
    reduction(+                                                                                   \
              : penergy)
  {

#ifdef _OPENMP
    offset = -1; // indicator
#else
    offset = 0; // zeroed
#endif

#pragma omp for
    for (i = 0; i < n_particles; i++)
    {

#ifdef _OPENMP
      // accumulate the offset now that we know i
      if (offset == -1)
      {
        offset = 0;
        for (j = 0; j < i; j++)
        {
          offset += nlist->nlist_count[j];
        }
      }
#endif
      // iterate through neighbor list
      //  for (n = offset; n - offset < nlist->nlist_count[i]; n++)
      //  {
      //    j = nlist->nlist[n];
      for (j = i + 1; j < n_particles + params->n_ghost_particles; j++)
      {
        r = 0;
        // get distance vector
        // put in 0 position
        for (k = 0; k < n_dims; k++)
        {
          // diff = min_image_dist(positions[j * n_dims + k] - positions[i * n_dims + k], box_size[k]);
          diff = positions[j * n_dims + k] - positions[i * n_dims + k];
          r += diff * diff;
          force_vector[k] = diff;
        }

        if (r > rcut * rcut)
          continue;
        r = sqrt(r);
        // LJ force and potential
        force = lj(r, epsilon, sigma, &e) - lj_shift;
#ifdef DEBUG
        printf("F(%d <-> %d, %g) = %g\n", i, j, r, force);
#endif // DEBUG

#pragma omp critical(update_forces)
        for (k = 0; k < n_dims; k++)
        {
          forces[i * n_dims + k] += force / r * force_vector[k];
          if (j < params->n_particles)
            forces[j * n_dims + k] -= force / r * force_vector[k];
        }

        penergy += e - e_shift;
      }
      // TODO: Turn back on when nlsit is working
      // offset += nlist->nlist_count[i];
    }
  }

  return penergy;
}

void lj_free_forces(force_t *force)
{
  lj_parameters_t *lj_parameters = (lj_parameters_t *)force->parameters;
  nlist_parameters_t *nlist = lj_parameters->nlist;
  free_nlist(nlist);
  free(lj_parameters);
  free(force);
}

force_t *build_lj(double epsilon, double sigma, nlist_parameters_t *nlist)
{
  lj_parameters_t init = {.epsilon = epsilon, .sigma = sigma};
  lj_parameters_t *parameters = (lj_parameters_t *)malloc(sizeof(lj_parameters_t));
  memcpy(parameters, &init, sizeof(init));
  parameters->nlist = nlist;
  force_t *lj = (force_t *)malloc(sizeof(force_t));
  lj->gather = &lj_gather_forces;
  lj->free = &lj_free_forces;
  lj->parameters = parameters;
  return lj;
}
