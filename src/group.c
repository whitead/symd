#include "group.h"
#include "box.h"
#include <string.h>
#include <math.h>

void *action(SCALAR *g, SCALAR *output, SCALAR *data, unsigned int n_dims, SCALAR s)
{
    unsigned int i, j;
    memset(output, 0, sizeof(SCALAR) * n_dims);
    for (i = 0; i < n_dims; i++)
    {
        for (j = 0; j < n_dims; j++)
            output[i] += data[j] * g[i * (n_dims + 1) + j];
        // w coord
        output[i] += s * g[i * (n_dims + 1) + j];
    }
}

void *fold_particles(run_params_t *params, SCALAR *positions)
{
    group_t *group = params->box->group;
    unsigned int n_dims = params->n_dims;
    const unsigned int p = params->n_particles;
    unsigned int i, j, k, l;
    // update scaled and wrap
    for (i = 0; i < p; i++)
        scale_wrap_coords(&params->scaled_positions[i * n_dims], &positions[i * n_dims], params->box);

    // unfold and update
    for (i = 0; i < p; i++)
    {
        for (j = 0; j < group->size; j++)
        {
            // unfold scaled, store temporarily in positions

            action(group->members[j].g, &positions[j * p * n_dims + i * n_dims],
                   &params->scaled_positions[i * n_dims], n_dims, 1.0);

            // tile and unscale
            for (k = 0; k < params->box->n_tilings; k++)
            {
                for (l = 0; l < n_dims; l++)
                    positions[n_dims * (p * ((k + 1) * group->size + j) + i) + l] =
                        positions[j * p * n_dims + i * n_dims] + params->box->tilings[k * n_dims + l];
                unscale_coords(&positions[n_dims * (p * ((k + 1) * group->size + j) + i)],
                               &positions[n_dims * (p * ((k + 1) * group->size + j) + i)], params->box);
            }
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
