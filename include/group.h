#include <stdbool.h>
#include "params.h"

#ifndef GROUP_H_
#define GROUP_H_

// just in case we want to add
// more info, kept as a struct
typedef struct
{
    SCALAR g[(N_DIMS + 1) * (N_DIMS + 1)];
} g_t;

struct group_t
{
    const char *name;
    unsigned int size;
    unsigned int total_size;
    unsigned int n_gparticles;
    unsigned int dof;
    g_t *members;
    SCALAR projector[N_DIMS * N_DIMS * N_DIMS * N_DIMS];
    SCALAR ainv[N_DIMS * N_DIMS]; // inverse for lagrange constraints
    group_t *next;
};

void fold_particles(run_params_t *params, SCALAR *positions);
void fold_velocities(run_params_t *params, SCALAR *velocities);
void update_group(group_t *group, box_t *box);
void apply_constraints(run_params_t *params, SCALAR *positions, SCALAR *velocities);
#endif // GROUP_H_