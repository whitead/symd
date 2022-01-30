#include "group.h"
#include "min_image.h"

static inline double *action(double *g, double *output, double *data, unsigned int n_dims)
{
    unsigned int i, j;
    for (i = 0; i < n_dims; i++)
        for (j = 0; j < n_dims; j++)
            output[i] += data[j] * g[i * n_dims + j];
    return output;
}

void *fold_particles(run_params_t *params, double *particles, bool reduce)
{
    group_t *group = params->group;
    unsigned int n_dims = params->n_dims, n_particles = params->n_particles + params->n_ghost_particles;
    const unsigned int p = n_particles / group->size;
    unsigned int i, j, k;
    for (i = 0; i < p; i++)
    {
        // fold to common partition
        if (reduce)
        {
            for (j = 0; j < group->size; j++)
            {
                // need to be wrapped
                // non-ghost should be wrapped in integrate
                for (k = 0; k < n_dims; k++)
                    particles[j * p * n_dims + i * n_dims + k] = wrap(particles[j * p * n_dims + i * n_dims + k], params->box_size[k]);
                action(group->members[j].i, &particles[i * n_dims], &particles[j * p * n_dims + i * n_dims], n_dims);
            }
            // average
            for (k = 0; k < n_dims; k++)
                particles[i * n_dims + k] /= group->size;
        }
        // unfold
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