#include "force.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

typedef struct
{
  const unsigned int n_types; 
  const double *epsilon;
  const double *sigma;
} lj_parameters_t;

typedef struct
{
  const double epsilon;
  const double sigma;
  nlist_parameters_t *nlist;
} nlj_parameters_t;

static inline void lj(int n_types, double r, double epsilon, double sigma, double result1[2*n_types], double result2[2*n_types])
{
  double ri6[n_types*n_types] = pow(sigma / r, 6);
  double ri12[n_types*n_types] = ri6 * ri6;
  result1 = 4 * epsilon * (6 * ri6 * sigma / r - 12 * ri12 * sigma / r);
  result2 = 4 * epsilon * (ri12 - ri6);
}

double nlj_gather_forces(int n_types, run_params_t *params, double *positions, double *forces)
{
  unsigned int n_dims = N_DIMS;
  unsigned int n_particles = params->n_particles;
  double *box_size = params->box->box_size;
  force_t *force_p = params->force_parameters;
  nlj_parameters_t *parameters = (nlj_parameters_t *)force_p->parameters;
  const double epsilon = parameters->epsilon;
  const double sigma = parameters->sigma;
  nlist_parameters_t *nlist = parameters->nlist;

  // update neighbor list
  // TODO: Turn back on
  update_nlist(positions, box_size, n_dims, n_particles, params->n_ghost_particles, nlist);

  unsigned int i, j, k, n;
  int offset;
  double penergy = 0;
  double r, force, diff, e, lj_result1[2*n_types], lj_result2[2*n_types];
  double force_vector[n_dims];
  double rcut = sqrt(nlist->rcut);
  double e_shift, lj_shift;
  lj(n_types, rcut, epsilon, sigma, lj_result1, lj_result2);
  lj_shift = lj_result1;
  e_shift = lj_result2;

#ifdef DEBUG
  check_nlist(params, nlist, positions, rcut);
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
          // TODO: turn back on
          offset += nlist->nlist_count[j];
        }
      }
#endif
      // iterate through neighbor list
      for (n = offset; n - offset < nlist->nlist_count[i]; n++)
      {
        j = nlist->nlist[n];
        r = 0;
        // get distance vector
        // put in 0 position
        for (k = 0; k < n_dims; k++)
        {
          diff = positions[j * n_dims + k] - positions[i * n_dims + k];
          r += diff * diff;
          force_vector[k] = diff;
        }

        if (r > rcut * rcut)
          continue;
        r = sqrt(r);
        // LJ force and potential
        lj(n_types, r, epsilon, sigma, lj_result1, lj_result2);
        force = lj_result1 - lj_shift;
        e = lj_result2 - e_shift;

#pragma omp critical(update_forces)
        for (k = 0; k < n_dims; k++)
        {
          forces[i * n_dims + k] += force / r * force_vector[k];
          if (j < params->n_particles)
            forces[j * n_dims + k] -= force / r * force_vector[k];
        }

        penergy += e;
      }
      offset += nlist->nlist_count[i];
    }
  }

  return penergy;
}

double lj_gather_forces(run_params_t *params, double *positions, double *forces)
{
  unsigned int n_dims = N_DIMS;
  unsigned int n_particles = params->n_particles;
  //unsigned int *types = params->types;
  force_t *force_p = params->force_parameters;
  lj_parameters_t *parameters = (lj_parameters_t *)force_p->parameters;
  const double *epsilon = parameters->epsilon;
  const double *sigma = parameters->sigma;

  unsigned int i, j, k;
  double penergy = 0;
  double r, force, diff, e, lj_result1, lj_result2;
  double force_vector[n_dims];
  // TODO better
  double *rcut = sigma * 3.5;
  double *e_shift, *lj_shift;
  lj(n_types, rcut, epsilon, sigma, lj_result1, lj_result2)
  lj_shift = lj_result1;
  e_shift = lj_result2;

// zero forces
#pragma omp parallel for
  for (i = 0; i < n_particles; i++)
    for (k = 0; k < n_dims; k++)
      forces[i * n_dims + k] = 0;

      // iterate through all particles

#pragma omp parallel default(shared) private(i, j, k, r, force, e, force_vector, diff) \
    reduction(+                                                                        \
              : penergy)
  {

#pragma omp for
    for (i = 0; i < n_particles; i++)
    {

      // iterate through remaining particles list
      for (j = i + 1; j < n_particles + params->n_ghost_particles; j++)
      {
        r = 0;
        // get distance vector
        // put in 0 position
        for (k = 0; k < n_dims; k++)
        {
          diff = positions[j * n_dims + k] - positions[i * n_dims + k];
          r += diff * diff;
          force_vector[k] = diff;
        }

        if (r > rcut * rcut)
          continue;
        r = sqrt(r);
        // LJ force and potential
        lj(n_types, r, epsilon, sigma, lj_result1, lj_result2)
        force = lj_result1 - lj_shift;
        e = lj_result2 - e_shift;

#pragma omp critical(update_forces)
        for (k = 0; k < n_dims; k++)
        {
          forces[i * n_dims + k] += force / r * force_vector[k];
          if (j < params->n_particles)
            forces[j * n_dims + k] -= force / r * force_vector[k];
        }
        penergy += e / 2;
        if (j < params->n_particles)
          penergy += e / 2;
      }
    }
  }

  return penergy;
}

void nlj_free_forces(force_t *force)
{
  nlj_parameters_t *lj_parameters = (nlj_parameters_t *)force->parameters;
  nlist_parameters_t *nlist = lj_parameters->nlist;
  free_nlist(nlist);
  free(lj_parameters);
  free(force);
}

force_t *build_nlj_types(int n_types, double epsilon, double sigma, nlist_parameters_t *nlist)
{
  double *sigma_new = (lj_parameters_t *)malloc(sizeof(n_types*n_types));
    double *epsilon_new = (lj_parameters_t *)malloc(sizeof(n_types*n_types));
    int index;
    for (int i=0; i<n_types; i++)
    {
        for (int j=0; j<n_types; j++)
        {
            index = i*n_types + j;
            sigma_new[index] = (sigma[i]*sigma[j])/2;
            epsilon_new[index] = (epsilon[i]*epsilon[j])/2
        }
    }
  nlj_parameters_t init = {.epsilon = epsilon_new, .sigma = sigma_new};
  nlj_parameters_t *parameters = (nlj_parameters_t *)malloc(sizeof(nlj_parameters_t));
  memcpy(parameters, &init, sizeof(init));
  parameters->nlist = nlist;
  force_t *lj = (force_t *)malloc(sizeof(force_t));
  lj->gather = &nlj_gather_forces;
  lj->free = &nlj_free_forces;
  lj->parameters = parameters;
  return lj;

void lj_free_forces(force_t *force)
{
  lj_parameters_t *lj_parameters = (lj_parameters_t *)force->parameters;
  free(lj_parameters);
  free(force);
}

force_t *build_lj_types(int n_types, double epsilon, double sigma)
{
    double *sigma_new = (lj_parameters_t *)malloc(sizeof(n_types*n_types));
    double *epsilon_new = (lj_parameters_t *)malloc(sizeof(n_types*n_types));
    int index;
    for (int i=0; i<n_types; i++)
    {
        for (int j=0; j<n_types; j++)
        {
            index = i*n_types + j;
            sigma_new[index] = (sigma[i]*sigma[j])/2;
            epsilon_new[index] = (epsilon[i]*epsilon[j])/2
        }
    }
    lj_parmeters_t init = {.epsilon = epsilon_new, .sigma = sigma_new};
    lj_parameters_t *parameters = (lj_parameters_t *)malloc(sizeof(lj_parameters_t));
    memcpy(parameters, &init, sizeof(init));
    force_t *lj = (force_t *)malloc(sizeof(force_t));
    lj->gather = &lj_gather_forces;
    lj->free = &lj_free_forces;
    lj->parameters = parameters;
    return lj;
}