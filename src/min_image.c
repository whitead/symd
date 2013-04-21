#include "min_image.h"
#include <stdio.h>

double round(double number)
{
  return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

double min_image_dist(double dx, double img) {
  return dx - round(dx / img) *img;
}

double wrap(double x, double img) {
  return x - floor(x /  img) * img;
}
