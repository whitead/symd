#include "params.h"

#ifndef MIN_IMAGE_H_
#define MIN_IMAGE_H_

/*
 * Minimum image distance functions
 *
 */
double round(double number);

double min_image_dist(double dx, double img);

double wrap(double x, double img);

double volume(double *box, unsigned int n_dims);

int try_rescale(run_params_t *params, double *positions, double *penergy, double *forces);

void tile(run_params_t *params, double *positions);

struct box_t {
    double* box_size;
    group_t* group;
};

#endif
