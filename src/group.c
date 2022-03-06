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
        // output[i] = fmod(output[i], 1.0);// ????
    }
}

static unsigned int _fold_particles(run_params_t *params, group_t *group, SCALAR *positions, unsigned int i_offset, unsigned int index_offset)
{
    const unsigned int p = params->n_particles;
    const unsigned int nc = params->n_cell_particles;
    SCALAR temp[N_DIMS];
    unsigned int i, j, k, l, index = index_offset;

    // unfold and update
#pragma omp parallel for default(shared) private(i, j, k, l, temp)
    for (j = 0; j < group->size; j++)
    {
        for (i = 0; i < group->n_gparticles; i++)
        {
            // unfold scaled, store temporarily in positions
            action(group->members[j].g, &positions[index * N_DIMS],
                   &params->scaled_positions[(i + i_offset) * N_DIMS], N_DIMS, 1.0);
            // tile and unscale
            for (k = 0; k < params->box->n_tilings; k++)
            {
                for (l = 0; l < N_DIMS; l++)
                    temp[l] = positions[index * N_DIMS + l] + params->box->tilings[k * N_DIMS + l];
                unscale_coords(&positions[nc * N_DIMS * k + index * N_DIMS],
                               temp, params->box);
            }
            // unscale the non-tiled folded scaled coordinates
            unscale_coords(temp,
                           &positions[index * N_DIMS], params->box);
            memcpy(&positions[index * N_DIMS], temp, N_DIMS * sizeof(SCALAR));
            index++;
        }
    }
    return index;
}

void fold_particles(run_params_t *params, SCALAR *positions)
{
    unsigned int i, j;
    // update scaled and wrap
#pragma omp parallel for default(shared) private(i)
    for (i = 0; i < params->n_particles; i++)
        scale_wrap_coords(&params->scaled_positions[i * N_DIMS], &positions[i * N_DIMS], params->box);
    // i is real particle, j is folded particle
    i = 0, j = 0;
    for (group_t *group = params->box->group; group != NULL; group = group->next)
    {
        j += _fold_particles(params, group, positions, i, j);
        i += group->n_gparticles;
    }
}

void fold_velocities(run_params_t *params, SCALAR *velocities)
{
    group_t *group = params->box->group;
    const unsigned int p = params->n_particles;
    SCALAR temp[N_DIMS];
    unsigned int i, j, k, l;

    // zero other velocities
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
    // TODO: Figure out why this segfaults someday.
    if (g == NULL)
        return;
    free(g->members);
    free_group(g->next);
    free(g);
}
