#include "group.h"
#include "box.h"
#include <string.h>
#include <math.h>

void action(SCALAR *g, SCALAR *output, SCALAR *data, unsigned int n_dims, SCALAR s)
{
    unsigned int i, j;
    memset(output, 0, sizeof(SCALAR) * n_dims);
    for (i = 0; i < n_dims; i++)
    {
        for (j = 0; j < n_dims; j++)
            output[i] += data[j] * g[i * (n_dims + 1) + j];
        // w coord
        output[i] += s * g[i * (n_dims + 1) + j];
        // output[i] = fmod(output[i], 1.0); ????
    }
}

void fold_particles(run_params_t *params, SCALAR *positions)
{
    group_t *group = params->box->group;
    const unsigned int p = params->n_particles;
    SCALAR temp[N_DIMS];
    unsigned int i, j, k, l;
// update scaled and wrap
#pragma omp parallel for default(shared) private(i)
    for (i = 0; i < p; i++)
        scale_wrap_coords(&params->scaled_positions[i * N_DIMS], &positions[i * N_DIMS], params->box);

        // unfold and update
#pragma omp parallel for default(shared) private(i, j, k, l, temp)
    for (i = 0; i < p; i++)
    {
        for (j = 0; j < group->size; j++)
        {
            // unfold scaled, store temporarily in positions

            action(group->members[j].g, &positions[j * p * N_DIMS + i * N_DIMS],
                   &params->scaled_positions[i * N_DIMS], N_DIMS, 1.0);
            // tile and unscale
            for (k = 0; k < params->box->n_tilings; k++)
            {
                for (l = 0; l < N_DIMS; l++)
                    temp[l] = positions[j * p * N_DIMS + i * N_DIMS + l] + params->box->tilings[k * N_DIMS + l];
                unscale_coords(&positions[N_DIMS * (p * ((k + 1) * group->size + j) + i)],
                               temp, params->box);
            }
            // unscale default non-tiling
            unscale_coords(temp,
                           &positions[j * p * N_DIMS + i * N_DIMS], params->box);
            memcpy(&positions[j * p * N_DIMS + i * N_DIMS], temp, N_DIMS * sizeof(SCALAR));
        }
    }
}

void fold_velocities(run_params_t *params, SCALAR *velocities)
{
    group_t *group = params->box->group;
    const unsigned int p = params->n_particles;
    SCALAR temp[N_DIMS];
    unsigned int i, j, k, l;

    //zero other velocities
    memset(&velocities[p], 0, params->n_ghost_particles * sizeof(SCALAR) * N_DIMS);

// unfold and update
#pragma omp parallel for default(shared) private(i, j, k, l, temp)
    for (i = 0; i < p; i++)
    {
        for (j = 1; j < group->size; j++)
        {
            action(group->members[j].g, &velocities[j * p * N_DIMS + i * N_DIMS],
                    &velocities[i * N_DIMS], N_DIMS, 0.0);
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
