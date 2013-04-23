#include "min_image.h"
#include <stdio.h>

double round(double number)
{
  return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

double min_image_dist(double dx, double img) {
#ifdef NO_PBC
  return dx;
#else
  return dx - round(dx / img) *img;
#endif
}

double wrap(double x, double img) {
#ifdef NO_PBC
  return x;
#else
  return x - floor(x /  img) * img;
#endif
}
