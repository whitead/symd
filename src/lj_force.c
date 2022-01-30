#include "force.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

static inline double lj(double r, double epsilon, double sigma)
{
  return 4 * epsilon * (6 * pow(sigma / r, 7) - 12 * pow(sigma / r, 13));
}

static inline double lj_trunc_shift(double r, double epsilon, double sigma, double rc, double shift)
{
  if (r >= rc)
  {
    return 0;
  }
  return lj(r, epsilon, sigma) - shift;
}

double lj_gather_forces(void *parameters, double *positions, double *forces, double *masses,
                        double *box_size, unsigned int n_dims, unsigned int n_particles)
{
  const double epsilon = ((Lj_Parameters *)parameters)->epsilon;
  const double sigma = ((Lj_Parameters *)parameters)->sigma;
  Nlist_Parameters *nlist = ((Lj_Parameters *)parameters)->nlist;

  //update neighbor list
  update_nlist(positions, box_size, n_dims, n_particles, nlist);

  unsigned int i, j, k, n;
  int offset;
  double penergy = 0;
  double r, force, diff;
  double force_vector[n_dims];
  double rcut = sqrt(nlist->rcut);
  double lj_shift = lj(rcut, epsilon, sigma);

//zero forces
#pragma omp parallel for
  for (i = 0; i < n_particles; i++)
    for (k = 0; k < n_dims; k++)
      forces[i * n_dims + k] = 0;

      //iterate through all particles

      //This seems strange at first,
      //but it really just distributes
      //the work across the threads
      //The only trick is that the neighborlists
      //are not easy to spread out across the workers,
      //hence the conditionals.

#pragma omp parallel default(shared) private(offset, n, i, j, k, r, force, force_vector, diff) \
    reduction(+                                                                                \
              : penergy)
  {

#ifdef _OPENMP
    offset = -1; //indicator
#else
    offset = 0; //zeroed
#endif

#pragma omp for
    for (i = 0; i < n_particles; i++)
    {

#ifdef _OPENMP
      //accumulate the offset now that we know i
      if (offset == -1)
      {
        offset = 0;
        for (j = 0; j < i; j++)
        {
          offset += nlist->nlist_count[j];
        }
      }
#endif

      //iterate through neighbor list
      for (n = offset; n - offset < nlist->nlist_count[i]; n++)
      {
        j = nlist->nlist[n];
        r = 0;

        //distance between particles
        for (k = 0; k < n_dims; k++)
        {
          diff = min_image_dist(positions[j * n_dims + k] - positions[i * n_dims + k], box_size[k]);
          r += diff * diff;
          force_vector[k] = diff;
        }

        r = sqrt(r);
        //LJ force and potential
        force = lj_trunc_shift(r, epsilon, sigma, rcut, lj_shift);

#ifdef DEBUG
        printf("F(%d - %d, %g) = %g\n", i, j, r, force);
#endif //DEBUG

#pragma omp critical(update_forces)
        for (k = 0; k < n_dims; k++)
        {
          forces[i * n_dims + k] += force / r * force_vector[k];
          forces[j * n_dims + k] -= force / r * force_vector[k];
        }

        penergy += 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6)) - 4 * epsilon * (pow(sigma / rcut, 12) - pow(sigma / rcut, 6));
      }

      offset += nlist->nlist_count[i];
    }
  }

  return (penergy);
}

void lj_free_forces(void *parameters)
{

  Lj_Parameters *lj_parameters = (Lj_Parameters *)parameters;
  Nlist_Parameters *nlist = ((Lj_Parameters *)lj_parameters)->nlist;
  free_nlist(nlist);
  free(lj_parameters);
}

force_t *build_lj(double epsilon, double sigma, Nlist_Parameters *nlist)
{
  Lj_Parameters init = {.epsilon = epsilon, .sigma = sigma};
  Lj_Parameters *parameters = (Lj_Parameters *)malloc(sizeof(Lj_Parameters));
  memcpy(parameters, &init, sizeof(init));
  parameters->nlist = nlist;
  force_t *lj = (force_t *)malloc(sizeof(force_t));
  lj->gather = &lj_gather_forces;
  lj->free = &lj_free_forces;
  lj->parameters = parameters;
  return lj;
}
