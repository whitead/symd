#include <math.h>

#ifndef MIN_IMAGE_H_
#define MIN_IMAGE_H_

/*
 * Minimum image distance functions
 *
 */
double round(double number);

double min_image_dist(double dx, double img);

double wrap(double x, double img);

#endif
