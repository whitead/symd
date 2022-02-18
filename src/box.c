#include "box.h"
#include "force.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>

static int sign(char x)
{
  return (x > 0) - (x < 0);
}

static int levi_civta(unsigned char index[N_DIMS])
{
  int p = 1;
  unsigned int i, j;
  for (i = 0; i < N_DIMS; i++)
    for (j = i + 1; j < N_DIMS; j++)
      p *= sign((char)(index[j]) - (char)(index[i]));
  return p;
}

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

// recursively find it. i refers to which vector
static SCALAR rvolume(SCALAR *b_vectors, SCALAR v, unsigned int i, unsigned char index[N_DIMS])
{
  if (i == N_DIMS)
    return levi_civta(index) * v;
  SCALAR vi = 0;
  unsigned int j;
  for (j = 0; j < N_DIMS; j++)
  {
    index[i] = j;
    vi += rvolume(b_vectors, b_vectors[j * N_DIMS + i] * v, i + 1, index);
  }
  return vi;
}

SCALAR volume(box_t *box)
{
  // https://math.stackexchange.com/a/1606559 -- Not sure if I believe it?
  unsigned char index[N_DIMS];
  return rvolume(box->b_vectors, 1, 0, index);
}

int try_rescale(run_params_t *params, SCALAR *positions, SCALAR *penergy, SCALAR *forces)
{
  unsigned int i, j;
  SCALAR newV, oldV = volume(params->box);
  SCALAR unorm_b_vectors[N_DIMS * N_DIMS];
  box_t *new_box = NULL;
  memcpy(unorm_b_vectors, params->box->b_vectors, sizeof(SCALAR) * N_DIMS * N_DIMS);

  // make random step along some sides
  // scale by .1% at most because it sounds reasonable
  char success = 0;
  while (!success)
  {
    for (i = 0; i < N_DIMS * N_DIMS; i++)
    {
      if (gsl_rng_uniform(params->rng) < 1.0 / n_dims)
      {
        // between 99% and 101%
        unorm_b_vectors[i] *= 1.0 + gsl_rng_uniform(params->rng) * 0.02 - 0.01;
        success = 1;
      }
    }
  }

#ifdef DEBUG
  printf("old box: ");
  for (i = 0; i < N_DIMS * N_DIMS; i++)
  {
    printf("%g ", params->box->b_vectors[i]);
  }
  printf("\n");
  printf("Proposed new box: ");
  for (i = 0; i < N_DIMS * N_DIMS; i++)
  {
    printf("%g ", unorm_b_vectors[i]);
  }
  printf("\n");
#endif

  new_box = make_box(unorm_b_vectors, params->box->group, N_DIMS, (params->box->n_tilings - 1) / 2);
  newV = volume(new_box, n_dims);

  // rescale coordinates
  // go to scaled
  for (i = 0; i < params->n_particles; i++)
    scale_wrap_coords(&params->scaled_positions[i * N_DIMS], &positions[i * N_DIMS], params->box);
  // unscale with new box
  for (i = 0; i < params->n_particles; i++)
    unscale_coords(&positions[i * N_DIMS], &params->scaled_positions[i * N_DIMS], new_box);

  // now fold them
  fold_particles(params, positions);

  // get energy
  SCALAR new_energy = params->force_parameters->gather(params, positions, forces);

  // accept/reject
  SCALAR mhc = (new_energy - *penergy) + params->pressure * (newV - oldV) -
               (params->n_particles + 1) * params->temperature * log(newV / oldV);

  if (gsl_rng_uniform(params->rng) < exp(-mhc / params->temperature))
  {
#ifdef DEBUG
    printf("Accepted with MHC %g\n", mhc);
#endif
    // accepted
    // TODO: rebuild cells in nlist
    *penergy = new_energy;
    params->box->group = NULL;
    free_box(params->box);
    params->box = new_box;
    return 1;
  }
  else
  {
#ifdef DEBUG
    printf("Rejected with MHC %g\n", mhc);
#endif
    // undo rescale coordinates
    // go to scaled
    for (i = 0; i < params->n_particles; i++)
      scale_wrap_coords(&params->scaled_positions[i * N_DIMS], &positions[i * N_DIMS], new_box);
    // unscale with old box
    for (i = 0; i < params->n_particles; i++)
      unscale_coords(&positions[i * N_DIMS], &params->scaled_positions[i * N_DIMS], params->box);

    // now fold them
    fold_particles(params, positions);
    new_box->group = NULL;
    free_box(new_box);
    return 0;
  }
}