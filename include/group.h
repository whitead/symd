#include <stdbool.h>
#include "params.h"

#ifndef GROUP_H_
#define GROUP_H_

// just in case we want to add
// more info, kept as a struct
typedef struct
{
    SCALAR *g;
} g_t;

struct group_t
{
    const char *name;
    unsigned int size;
    unsigned int total_size;
    unsigned int n_gparticles;
    unsigned int dof;
    g_t *members;
    SCALAR *projector;
    group_t *next;
};

void fold_particles(run_params_t *params, SCALAR *positions);
void fold_velocities(run_params_t *params, SCALAR *velocities);
void free_group(group_t *g);
#endif // GROUP_H_