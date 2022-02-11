#include <stdbool.h>
#include "params.h"

#ifndef GROUP_H_
#define GROUP_H_

typedef struct
{
    double *g;
    double *i;
    int tiling;
} g_t;

struct group_t
{
    const char *name;
    unsigned int size;
    unsigned int tiling_start;
    g_t *members;
    double* projector;
};

void *fold_particles(run_params_t *group, double *particles, double *velocities, bool reduce);
void free_group(group_t *g);
#endif // GROUP_H_