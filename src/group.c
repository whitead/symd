#include "group.h"
#include "box.h"
#include <math.h>

static inline double *action(double *g, double *output, double *data, double *box, unsigned int n_dims, double s)
{
    unsigned int i, j;
    for (i = 0; i < n_dims; i++)
    {
        for (j = 0; j < n_dims; j++)
        {
            // printf("output[i(%d)] (%f) = data[j] (%f) * g[i * (n_dims + 1) + j (%d)] (%f)\n", i, output[i], data[j], i * (n_dims + 1) + j, g[i * (n_dims + 1) + j]);
            output[i] += data[j] * g[i * (n_dims + 1) + j];
        }
        // w coord
        // printf("output[i(%d)] (%f) += s (%f) * g[i * (n_dims + 1) + j] (%f) * box[i] (%f);\n", i, output[i], s, g[i * (n_dims + 1) + j], box[i]);
        output[i] += s * g[i * (n_dims + 1) + j] * box[i];
    }
    return output;
}

void *fold_particles(run_params_t *params, double *particles, bool reduce, bool vector)
{
    group_t *group = params->group;
    unsigned int n_dims = params->n_dims;
    double r0, r;
    const unsigned int p = params->n_particles;
    unsigned int i, j, k;
    // TODO: OpenMP
    for (i = 0; i < p; i++)
    {
        // fold to common partition
        if (reduce)
        {
            for (j = 0; j < group->size; j++)
            {
                action(group->members[j].i,
                       &particles[i * n_dims], &particles[j * p * n_dims + i * n_dims],
                       params->box_size, n_dims, vector ? 0.0 : 1.0);
            }
            // average
            for (k = 0; k < n_dims; k++)
                particles[i * n_dims + k] /= group->size;
        }

        // get minimum image
        if (!vector)
        {
            // unfold
            for (j = group->tiling_start; j < group->size; j++)
            {
                memset(&particles[j * p * n_dims + i * n_dims], 0, sizeof(double) * n_dims);
                action(group->members[j].g, &particles[j * p * n_dims + i * n_dims],
                       &particles[i * n_dims], params->box_size, n_dims, vector ? 0.0 : 1.0);
                // printf("particles[j (%d) * p * n_dims + i * n_dims] = %f %f\n", j, particles[j * p * n_dims + i * n_dims], particles[j * p * n_dims + i * n_dims + 1]);
            }

            r0 = 0;
            for (k = 0; k < n_dims; k++)
                r0 += pow(particles[i * n_dims + k], 2);
            for (j = group->tiling_start; j < group->size; j++)
            {
                // check if minimum image
                r = 0;
                for (k = 0; k < n_dims; k++)
                    r += pow(particles[j * p * n_dims + i * n_dims + k], 2);
                if (r < r0)
                {
                    r0 = r;
                    memcpy(&particles[i * n_dims], &particles[j * p * n_dims + i * n_dims], sizeof(double) * n_dims);
                }
            }
        }

        // unfold
        for (j = 1; j < group->size; j++)
        {
            memset(&particles[j * p * n_dims + i * n_dims], 0, sizeof(double) * n_dims);
            action(group->members[j].g, &particles[j * p * n_dims + i * n_dims],
                   &particles[i * n_dims], params->box_size, n_dims, vector ? 0.0 : 1.0);
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