#include "box.h"
#include "force.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

static void
invert_matrix(SCALAR *data, SCALAR* inv_data)
{

    gsl_matrix_view matrix = gsl_matrix_view_array( data, N_DIMS, N_DIMS );
    gsl_permutation *p = gsl_permutation_alloc(N_DIMS);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(&matrix.matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix_view inv = gsl_matrix_view_array( inv_data, N_DIMS, N_DIMS );
    gsl_linalg_LU_invert(&matrix.matrix, p, &inv.matrix);

    gsl_permutation_free(p);
    
}


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
  unsigned int i, j;
  memset(dest, 0, N_DIMS * sizeof(SCALAR));
  for (i = 0; i < N_DIMS; i++)
  {
    for (j = 0; j < N_DIMS; j++)
      dest[i] += src[j] * box->ib_vectors[i * N_DIMS + j];
    dest[i] = fmod(dest[i], 1.0);
  }
}

void unscale_coords(SCALAR *dest, SCALAR *src, box_t *box)
{
  unsigned int i, j;
  memset(dest, 0, N_DIMS * sizeof(SCALAR));
  for (i = 0; i < N_DIMS; i++)
    for (j = 0; j < N_DIMS; j++)
      dest[i] += src[j] * box->b_vectors[i * N_DIMS + j];
}

void tiling(int *result, unsigned int cur_dim, unsigned int n)
{
  unsigned int i, j;
  for (i = 0; i < n * 2 + 1; i++)
    result[cur_dim * N_DIMS + i] = (int)(i)-n;
}

void free_box(box_t *box)
{
  if (box->group)
    free_group(box->group);
  if (box->unorm_b_vectors)
    free(box->unorm_b_vectors);
  free(box->box_size);
  free(box->b_vectors);
  free(box->ib_vectors);
  free(box->tilings);
  free(box);
}

box_t *make_box(SCALAR *unorm_b_vectors, group_t *group, unsigned int images)
{

  box_t *box = (box_t *)malloc(sizeof(box_t));
  box->box_size = (SCALAR *)calloc(N_DIMS, sizeof(SCALAR));
  box->b_vectors = (SCALAR *)calloc(N_DIMS * N_DIMS, sizeof(SCALAR));
  box->ib_vectors = (SCALAR *)malloc(N_DIMS * N_DIMS * sizeof(SCALAR));
  box->unorm_b_vectors = unorm_b_vectors;
  box->group = group;
  box->n_images = images;

  unsigned int i, j, k, l;

  // project to get valid basis vector
  if (group)
    for (i = 0; i < N_DIMS * N_DIMS; i++)
      for (j = 0; j < N_DIMS * N_DIMS; j++)
        box->b_vectors[i] += unorm_b_vectors[j] * group->projector[i * N_DIMS * N_DIMS + j];
  else
    memcpy(box->b_vectors, unorm_b_vectors, sizeof(SCALAR) * N_DIMS * N_DIMS);

  // compute their inverse
  invert_matrix(box->b_vectors, box->ib_vectors);

  // get max box size in each dimension
  for (i = 0; i < N_DIMS; i++)
    for (j = 0; j < N_DIMS; j++)
      box->box_size[i] = fmax(box->box_size[i], box->b_vectors[i * N_DIMS + j]);

#ifdef DEBUG
  printf("Box is size ");
  for (i = 0; i < N_DIMS; i++)
    printf("%f ", box->box_size[i]);
  printf("\nUnormed vectors (column wise):\n");
  for (i = 0; i < N_DIMS; i++)
  {
    for (j = 0; j < N_DIMS; j++)
      printf("%f ", unorm_b_vectors[i * N_DIMS + j]);
    printf("\n");
  }
  printf("Projected vectors:\n");
  for (i = 0; i < N_DIMS; i++)
  {
    for (j = 0; j < N_DIMS; j++)
      printf("%f ", box->b_vectors[i * N_DIMS + j]);
    printf("\n");
  }
  printf("Inverse Projected vectors:\n");
  for (i = 0; i < N_DIMS; i++)
  {
    for (j = 0; j < N_DIMS; j++)
      printf("%f ", box->ib_vectors[i * N_DIMS + j]);
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
  unsigned int ntot = pow(n, N_DIMS);
  // used to omit origin
  unsigned int past_zero = 0;
  int v;
  box->tilings = malloc(sizeof(int) * (ntot - 1) * N_DIMS);
  for (i = 0; i < ntot; i++)
  {
    l = i;
    if (i == ntot / 2)
    {
      past_zero = 1;
      continue;
    }
    for (j = 0; j < N_DIMS; j++)
    {
      v = l / pow(n, N_DIMS - j - 1);
      l -= v * pow(n, N_DIMS - j - 1);
      box->tilings[(i - past_zero) * N_DIMS + j] = v - images;
    }
  }
  box->n_tilings = ntot - 1;
#ifdef DEBUG
  printf("tilings:\n");
  for (i = 0; i < box->n_tilings; i++)
  {
    printf("%d: ", i);
    for (j = 0; j < N_DIMS; j++)
      printf("%d ", box->tilings[i * N_DIMS + j]);
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
  SCALAR *unorm_b_vectors = (SCALAR *)malloc(sizeof(SCALAR) * N_DIMS * N_DIMS);
  box_t *new_box = NULL;
  memcpy(unorm_b_vectors, params->box->unorm_b_vectors, sizeof(SCALAR) * N_DIMS * N_DIMS);

  // make random step along some sides
  // scale by .1% at most because it sounds reasonable
  char success = 0;
  while (!success)
  {
    for (i = 0; i < N_DIMS * N_DIMS; i++)
    {
      if (unorm_b_vectors[i] && gsl_rng_uniform(params->rng) < 1.0 / N_DIMS)
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

  new_box = make_box(unorm_b_vectors, params->box->group, params->box->n_images);
  newV = volume(new_box);

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
    printf("Accepted with MHC %g (delta E = %g, delta V = %g)\n", mhc, new_energy - *penergy, newV - oldV);
#endif // DEBUG
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
    printf("Rejected with MHC %g (delta E = %g, delta V = %g)\n", mhc, new_energy - *penergy, newV - oldV);
#endif // DEBUG
    // undo rescale coordinates
    // go to scaled
    for (i = 0; i < params->n_particles; i++)
      scale_wrap_coords(&params->scaled_positions[i * N_DIMS], &positions[i * N_DIMS], new_box);
    // unscale with old box
    for (i = 0; i < params->n_particles; i++)
      unscale_coords(&positions[i * N_DIMS], &params->scaled_positions[i * N_DIMS], params->box);

    // now fold them
    fold_particles(params, positions);

    // reset forces
    new_energy = params->force_parameters->gather(params, positions, forces);
    // just incase roundtrip rescale-rescale changes things
    *penergy = new_energy;
    new_box->group = NULL;
    free_box(new_box);
    return 0;
  }
}