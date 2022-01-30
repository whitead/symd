#include <stdbool.h>
#include "params.h"

#ifndef GROUP_H_
#define GROUP_H_

typedef struct
{
    double *g;
    double *i;
} g_t;

struct group_t
{
    const char *name;
    unsigned int size;
    g_t *members;
};

void *fold_particles(run_params_t *group, double *particles, bool reduce);
void free_group(group_t *g);
#endif // GROUP_H_