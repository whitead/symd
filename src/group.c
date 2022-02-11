#include "group.h"
#include "box.h"
#include <string.h>
#include <math.h>

double *action(double *g, double *output, double *data, double *basis_vectors, unsigned int n_dims, double s)
{
    unsigned int i, j;
    for (i = 0; i < n_dims; i++)
    {
        for (j = 0; j < n_dims; j++)
        {
            // printf("output[i(%d)] (%f) = data[j] (%f) * g[i * (n_dims + 1) + j (%d)] (%f)\n", i,
            //        output[i], data[j], i * (n_dims + 1) + j, g[i * (n_dims + 1) + j]);
            output[i] += data[j] * g[i * (n_dims + 1) + j];
        }
        // w coord
        // printf("output[i(%d)] (%f) += s (%f) * g[i * (n_dims + 1) + j] (%f) * box[i] (%f);\n",
        //        i, output[i], s, g[i * (n_dims + 1) + j], box[i]);
        output[i] += s * g[i * (n_dims + 1) + j] * box[i];
    }
    return output;
}

void *fold_particles(run_params_t *params, double *positions)
{
    group_t *group = params->box->group;
    unsigned int n_dims = params->n_dims;
    double r0, r;
    const unsigned int p = params->n_particles;
    unsigned int i, j, k;
    // update scaled and wrap
    for(i = 0; i < p; i++)
        scale_wrap_coords(&params->scaled_positions[i * n_dims], &positions[i * n_dims], params->box);

    // unfold and update
    for (i = 0; i < p; i++)
    {
        for (j = 0; j < group->size; j++)
        {
            //Update cartesian and unfold
            memset(&positions[j * p * n_dims + i * n_dims], 0, sizeof(double) * n_dims);
            action(group->members[j].g, &positions[j * p * n_dims + i * n_dims],
                   &positions[i * n_dims], params->box->box_size, n_dims, 1.0);
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
    free(g->projector);
    free(g);
}
