#include "box.h"
#include "force.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>

double round(double number)
{
  return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

double min_image_dist(double dx, double img)
{
  if (img)
    return dx - round(dx / img) * img;
  return dx;
}

double wrap(double x, double img)
{
  if (img)
  {
    double cx = x + img / 2;
    x = cx - floor(cx / img) * img;
    x = x - img / 2;
  }
  return x;
}

double volume(double *box, unsigned int n_dims)
{
  double v = 1.0;
  for (unsigned int i = 0; i < n_dims; i++)
    v *= box[i];
  return v;
}

int try_rescale(run_params_t *params, double *positions, double *penergy, double *forces)
{
  unsigned int i, j, n_dims = params->n_dims;
  double newV, oldV = volume(params->box->box_size, params->n_dims);
  double new_box[n_dims];
  memcpy(new_box, params->box->box_size, sizeof(double) * n_dims);

  // make random step along some sides
  if (params->box->kind == PBC_CUBIC || params->box->kind == GROUP_CUBIC)
  {
    double dx = 1.0 + gsl_rng_uniform(params->rng) * 0.02 - 0.01;
    for (i = 0; i < n_dims; i++)
      new_box[i] *= dx;
  }
  else if (params->box->kind == PBC)
  {
    // scale by .1% at most because it sounds reasonable
    int success = 0;
    while (!success)
    {
      for (i = 0; i < n_dims; i++)
      {
        if (gsl_rng_uniform(params->rng) < 1.0 / n_dims)
        {
          //between 99% and 101%
          new_box[i] *= 1.0 + gsl_rng_uniform(params->rng) * 0.02 - 0.01;
          success = 1;
        }
      }
    }
  }

#ifdef DEBUG
  printf("old box: ");
  for (i = 0; i < n_dims; i++)
  {
    printf("%g ", params->box->box_size[i]);
  }
  printf("\n");
  printf("Proposed new box: ");
  for (i = 0; i < n_dims; i++)
  {
    printf("%g ", new_box[i]);
  }
  printf("\n");
#endif

  newV = volume(new_box, n_dims);

  //rescale coordinates
  for (i = 0; i < params->n_ghost_particles + params->n_particles; i++)
    for (j = 0; j < n_dims; j++)
      positions[i * n_dims + j] *= new_box[j] / params->box->box_size[j];

  double new_energy = params->force_parameters->gather(params, positions, forces);

  // accept/reject
  double mhc = (new_energy - *penergy) + params->pressure * (newV - oldV) -
               (params->n_particles + 1) * params->temperature * log(newV / oldV);

  if (gsl_rng_uniform(params->rng) < exp(-mhc / params->temperature))
  {
#ifdef DEBUG
    printf("Accepted with MHC %g\n", mhc);
#endif
    //accepted
    //TODO: rebuild cells in nlist - shouldn't matter if cubic
    *penergy = new_energy;
    memcpy(params->box->box_size, new_box, sizeof(double) * n_dims);
    return 1;
  }
  else
  {
#ifdef DEBUG
    printf("Rejected with MHC %g\n", mhc);
#endif
    //undo rescale coordinates
    for (i = 0; i < params->n_ghost_particles + params->n_particles; i++)
      for (j = 0; j < n_dims; j++)
        positions[i * n_dims + j] /= new_box[j] / params->box->box_size[j];
    return 0;
  }
}