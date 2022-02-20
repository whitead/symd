#include "params.h"
#include "group.h"

#ifndef MIN_IMAGE_H_
#define MIN_IMAGE_H_

double min_image_dist(double dx, double img);

double wrap(double x, double img);

void scale_wrap_coords(SCALAR *dest, SCALAR *src, box_t *box);

void unscale_coords(SCALAR *dest, SCALAR *src, box_t *box);

double volume(box_t *box);

int try_rescale(run_params_t *params, double *positions, double *penergy, double *forces);

void tile(run_params_t *params, double *positions);

box_t *make_box(double *unorm_b_vectors, group_t *group, unsigned int n_tilings);
void free_box(box_t *box);

struct box_t
{
    double *cell_size;
    // for these -> columns are basis vectors!
    double *b_vectors;
    double *unorm_b_vectors;
    double *ib_vectors;
    int *tilings;
    unsigned int n_tilings;
    unsigned int n_images;
    group_t *group;
};

#endif
