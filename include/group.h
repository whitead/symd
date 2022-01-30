#ifndef GROUP_H_
#define GROUP_H_

typedef struct
{
    double *g;
    double *i;
} g_t;

typedef struct
{
    const char *name;
    unsigned int size;
    g_t *members;
} group_t;

void *fold_particles(group_t *group, double *particles, unsigned int n_dims, unsigned int n_particles);
void free_group(group_t *g);
#endif // GROUP_H_