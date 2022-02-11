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

double volume(double *box, unsigned int n_dims);

int try_rescale(run_params_t *params, double *positions, double *penergy, double *forces);

void tile(run_params_t *params, double *positions);

g_t *box_dist_start(box_t *box);
g_t *box_dist_next(box_t *box, g_t *g, double *dx, double *r);
unsigned int box_dist_size(box_t *box);

struct box_t
{
    double *box_size;
    double *b_vectors;
    double *unorm_b_vectors;
    unsigned int n_dims;
    group_t *group;
    box_e kind;
};

#endif
