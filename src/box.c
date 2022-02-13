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
  return dx - round(dx / img) * img;
}

double wrap(double x, double img)
{
  return x - floor(x / img) * img;
}

void scale_wrap_coords(SCALAR *dest, SCALAR *src, box_t *box)
{
  unsigned int n_dims = box->n_dims;
  unsigned int i, j;
  memset(dest, 0, n_dims * sizeof(SCALAR));
  for (i = 0; i < n_dims; i++)
  {
    for (j = 0; j < n_dims; j++)
      dest[i] += src[j] * box->ib_vectors[i * n_dims + j];
    dest[i] = fmod(dest[i], 1.0);
  }
}

void unscale_coords(SCALAR *dest, SCALAR *src, box_t *box)
{
  unsigned int n_dims = box->n_dims;
  unsigned int i, j;
  memset(dest, 0, n_dims * sizeof(SCALAR));
  for (i = 0; i < n_dims; i++)
    for (j = 0; j < n_dims; j++)
      dest[i] += src[j] * box->b_vectors[i * n_dims + j];
}

void tiling(int *result, unsigned int cur_dim, unsigned int n, unsigned n_dims)
{
  unsigned int i, j;
  for (i = 0; i < n * 2 + 1; i++)
    result[cur_dim * n_dims + i] = (int)(i)-n;
}

void free_box(box_t *box)
{
  if (box->group)
    free_group(box->group);
  free(box->box_size);
  free(box->b_vectors);
  free(box->ib_vectors);
  free(box->unorm_b_vectors);
  free(box->tilings);
  free(box);
}

box_t *make_box(SCALAR *unorm_b_vectors, group_t *group, unsigned int n_dims, unsigned int images)
{

  box_t *box = (box_t *)malloc(sizeof(box_t));
  box->box_size = (SCALAR *)calloc(n_dims, sizeof(SCALAR));
  box->b_vectors = (SCALAR *)calloc(n_dims * n_dims, sizeof(SCALAR));
  box->ib_vectors = (SCALAR *)malloc(n_dims * n_dims * sizeof(SCALAR));
  box->unorm_b_vectors = unorm_b_vectors;
  box->n_dims = n_dims;
  box->group = group;

  unsigned int i, j, k, l;

  // project to get valid basis vector
  if (group)
    for (i = 0; i < n_dims * n_dims; i++)
      for (j = 0; j < n_dims * n_dims; j++)
        box->b_vectors[i] += unorm_b_vectors[j] * group->projector[i * n_dims * n_dims + j];
  else
    box->b_vectors = unorm_b_vectors;

  // compute their inverse
  if (n_dims == 2)
  {
    double s = 1 / (box->b_vectors[0] * box->b_vectors[3] - box->b_vectors[1] * box->b_vectors[2]);
    box->ib_vectors[0] = box->b_vectors[3] * s;
    box->ib_vectors[1] = -box->b_vectors[1] * s;
    box->ib_vectors[2] = -box->b_vectors[2] * s;
    box->ib_vectors[3] = box->b_vectors[0] * s;
  }
  else
  {
    fprintf(stderr, "Finish implementing matrix invese\n");
    exit(1);
  }

  // get max box size in each dimension
  for (i = 0; i < n_dims; i++)
    for (j = 0; j < n_dims; j++)
      box->box_size[i] = fmax(box->box_size[i], box->b_vectors[i * n_dims + j]);

#ifdef DEBUG
  printf("Box is size ");
  for (i = 0; i < n_dims; i++)
    printf("%f ", box->box_size[i]);
  printf("\nUnormed vectors (column wise):\n");
  for (i = 0; i < n_dims; i++)
  {
    for (j = 0; j < n_dims; j++)
      printf("%f ", box->unorm_b_vectors[i * n_dims + j]);
    printf("\n");
  }
  printf("Projected vectors:\n");
  for (i = 0; i < n_dims; i++)
  {
    for (j = 0; j < n_dims; j++)
      printf("%f ", box->b_vectors[i * n_dims + j]);
    printf("\n");
  }
  printf("Inverse Projected vectors:\n");
  for (i = 0; i < n_dims; i++)
  {
    for (j = 0; j < n_dims; j++)
      printf("%f ", box->ib_vectors[i * n_dims + j]);
    printf("\n");
  }

#endif

  // build out tilings
  // goes through integers -> 0, 1, 2, ...
  // converts to base of images * 2 + 1
  // so if we want to having tilings at -2, -1, 0, 1, 2
  // we get that by using base 5 and an integer
  // like 50 would be 200 in base 5, which would be 0,-2,-2 for 3D
  unsigned int n = images * 2 + 1;
  unsigned int ntot = pow(n, n_dims);
  // used to omit origin
  unsigned int past_zero = 0;
  int v;
  box->tilings = malloc(sizeof(int) * (ntot - 1) * n_dims);
  for (i = 0; i < ntot; i++)
  {
    l = i;
    if (i == ntot / 2)
    {
      past_zero = 1;
      continue;
    }
    for (j = 0; j < n_dims; j++)
    {
      v = l / pow(n, n_dims - j - 1);
      l -= v * pow(n, n_dims - j - 1);
      box->tilings[(i - past_zero) * n_dims + j] = v - images;
    }
  }
  box->n_tilings = ntot - 1;
#ifdef DEBUG
  printf("tilings:\n");
  for (i = 0; i < box->n_tilings; i++)
  {
    printf("%d: ", i);
    for (j = 0; j < n_dims; j++)
      printf("%d ", box->tilings[i * n_dims + j]);
    printf("\n");
  }
#endif

  return box;
}

double volume(double *box, unsigned int n_dims)
{
  double v = 1.0;
  for (unsigned int i = 0; i < n_dims; i++)
    v *= box[i];
  return v;
}

// int try_rescale(run_params_t *params, double *positions, double *penergy, double *forces)
// {
//   unsigned int i, j, n_dims = N_DIMS;
//   double newV, oldV = volume(params->box->box_size, N_DIMS);
//   double new_box[n_dims];
//   memcpy(new_box, params->box->box_size, sizeof(double) * n_dims);

//   // make random step along some sides
//   if (params->box->kind == PBC_CUBIC || params->box->kind == GROUP_CUBIC)
//   {
//     double dx = 1.0 + gsl_rng_uniform(params->rng) * 0.02 - 0.01;
//     for (i = 0; i < n_dims; i++)
//       new_box[i] *= dx;
//   }
//   else if (params->box->kind == PBC)
//   {
//     // scale by .1% at most because it sounds reasonable
//     int success = 0;
//     while (!success)
//     {
//       for (i = 0; i < n_dims; i++)
//       {
//         if (gsl_rng_uniform(params->rng) < 1.0 / n_dims)
//         {
//           // between 99% and 101%
//           new_box[i] *= 1.0 + gsl_rng_uniform(params->rng) * 0.02 - 0.01;
//           success = 1;
//         }
//       }
//     }
//   }

// #ifdef DEBUG
//   printf("old box: ");
//   for (i = 0; i < n_dims; i++)
//   {
//     printf("%g ", params->box->box_size[i]);
//   }
//   printf("\n");
//   printf("Proposed new box: ");
//   for (i = 0; i < n_dims; i++)
//   {
//     printf("%g ", new_box[i]);
//   }
//   printf("\n");
// #endif

//   newV = volume(new_box, n_dims);

//   // rescale coordinates
//   for (i = 0; i < params->n_ghost_particles + params->n_particles; i++)
//     for (j = 0; j < n_dims; j++)
//       positions[i * n_dims + j] *= new_box[j] / params->box->box_size[j];

//   double new_energy = params->force_parameters->gather(params, positions, forces);

//   // accept/reject
//   double mhc = (new_energy - *penergy) + params->pressure * (newV - oldV) -
//                (params->n_particles + 1) * params->temperature * log(newV / oldV);

//   if (gsl_rng_uniform(params->rng) < exp(-mhc / params->temperature))
//   {
// #ifdef DEBUG
//     printf("Accepted with MHC %g\n", mhc);
// #endif
//     // accepted
//     // TODO: rebuild cells in nlist - shouldn't matter if cubic
//     *penergy = new_energy;
//     memcpy(params->box->box_size, new_box, sizeof(double) * n_dims);
//     return 1;
//   }
//   else
//   {
// #ifdef DEBUG
//     printf("Rejected with MHC %g\n", mhc);
// #endif
//     // undo rescale coordinates
//     for (i = 0; i < params->n_ghost_particles + params->n_particles; i++)
//       for (j = 0; j < n_dims; j++)
//         positions[i * n_dims + j] /= new_box[j] / params->box->box_size[j];
//     return 0;
//   }
// }