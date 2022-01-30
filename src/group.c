#include "group.h"
#include "util.h"

static inline double *action(double *g, double *output, double *data, unsigned int n_dims)
{
    unsigned int i, j;
    for (i = 0; i < n_dims; i++)
        for (j = 0; j < n_dims; j++)
            output[i] += data[j] * g[i * n_dims + j];
    return output;
}

void *fold_particles(group_t *group, double *particles, unsigned int n_dims, unsigned int n_particles)
{
    const unsigned int p = n_particles / group->size;
    unsigned int i, j;
    for (i = 0; i < p; i++)
    {
        for (j = 1; j < group->size; j++)
            action(group->members[j].i, &particles[i * n_dims], &particles[j * p * n_dims + i * n_dims], n_dims);
        for (j = 0; j < n_dims; j++)
            particles[i * n_dims + j] /= group->size;
        for (j = 1; j < group->size; j++)
        {
            memset(&particles[j * p * n_dims + i * n_dims], 0, sizeof(double) * n_dims);
            action(group->members[j].g, &particles[j * p * n_dims + i * n_dims], &particles[i * n_dims], n_dims);
        }
    }
}

void free_group(group_t *g)
{
    for (unsigned int i = 0; i < g->size; i++)
    {
        free(g->members[i].g);
        free(g->members[i].i);
    }
    free(g->members);
    free(g);
}