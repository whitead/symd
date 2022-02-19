#include <stdbool.h>
#include "params.h"

#ifndef GROUP_H_
#define GROUP_H_

typedef struct
{
    SCALAR *g;
    SCALAR *i;
} g_t;

struct group_t
{
    const char *name;
    unsigned int size;
    g_t *members;
    SCALAR *projector;
};

void fold_particles(run_params_t *params, SCALAR *positions);
void free_group(group_t *g);
#endif // GROUP_H_