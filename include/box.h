#include "params.h"
#include "group.h"

#ifndef MIN_IMAGE_H_
#define MIN_IMAGE_H_

typedef enum box_e_
{
    UNWRAPPED,
    PBC,
    PBC_CUBIC,
    GROUP,
    GROUP_CUBIC
} box_e;

/*
 * Minimum image distance functions
 *
 */
double
round(double number);

double min_image_dist(double dx, double img);

double wrap(double x, double img);

void scale_wrap_coords(SCALAR *dest, SCALAR *src, box_t *box);

void unscale_coords(SCALAR *dest, SCALAR *src, box_t *box);

double volume(double *box, unsigned int n_dims);

int try_rescale(run_params_t *params, double *positions, double *penergy, double *forces);

void tile(run_params_t *params, double *positions);

box_t *make_box(double *unorm_b_vectors, group_t *group, unsigned int n_dims, unsigned int n_tilings);

struct box_t
{
    double *box_size;
    double *b_vectors;
    double *ib_vectors;
    double *unorm_b_vectors;
    int *tilings;
    unsigned int n_dims;
    unsigned int n_tilings;
    group_t *group;
    box_e kind;
};

#endif
