#include "min_image.h"
#include <stdio.h>

double round(double number)
{
  return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

double min_image_dist(double dx, double img)
{
  if (img)
    return dx - round(dx / img) * img;
  return dx;
}

double wrap(double x, double img)
{
  if (img)
  {
    double cx = x + img / 2;
    x = cx - floor(cx / img) * img;
    return x - img / 2;
  }
  return x;
}
