#include "group.h"
#include "box.h"
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

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
        // output[i] = fmod(output[i], 1.0);
    }
}

static unsigned int _fold_particles(run_params_t *params, group_t *group, SCALAR *positions, unsigned int i_offset, unsigned int index_offset)
{
    const unsigned int p = params->n_particles;
    const unsigned int nc = params->n_cell_particles;
    SCALAR temp[N_DIMS];
    unsigned int i, j, k, l, index;

    // tile
    // asymmetric unit
    j = 0;
    index = i_offset;
#pragma omp parallel for default(shared) private(i, index, k, l, temp)
    for (i = 0; i < group->n_gparticles; i++, index++)
    {
        // tile and unscale
        for (k = 0; k < params->box->n_tilings; k++)
        {
            for (l = 0; l < N_DIMS; l++)
                temp[l] = params->scaled_positions[index * N_DIMS + l] + params->box->tilings[k * N_DIMS + l];
            unscale_coords(&positions[nc * N_DIMS * (k + 1) + index * N_DIMS],
                           temp, params->box);
        }
        // this is done to ensure wrapping of asymmetric particles
        unscale_coords(&positions[index * N_DIMS],
                       &params->scaled_positions[index * N_DIMS], params->box);
#ifdef DEBUG
        printf("Tiling only group %s: %d\n", group->name, i + i_offset);
        printf("us: %f %f -> s: %f %f\n", params->scaled_positions[index * N_DIMS],
               params->scaled_positions[index * N_DIMS + 1], positions[index * N_DIMS],
               positions[index * N_DIMS + 1]);
#endif
    }
    // unfold and tile
    // TODO: can we get it to use OpenMP?
    index = p + index_offset;
    for (j = 1; j < group->size; j++)
    {
        for (i = 0; i < group->n_gparticles; i++, index++)
        {
            // unfold scaled, store temporarily in positions
            action(group->members[j].g, &positions[index * N_DIMS],
                   &params->scaled_positions[(i + i_offset) * N_DIMS], N_DIMS, 1.0);
#ifdef DEBUG
            printf("Applying group %s:%d for %d -> %d\n", group->name, j, i + i_offset, index);
            printf("us: %f %f -> %f %f\n", params->scaled_positions[(i + i_offset) * N_DIMS],
                   params->scaled_positions[(i + i_offset) * N_DIMS + 1], positions[index * N_DIMS + 0],
                   positions[index * N_DIMS + 1]);
#endif
            // tile and unscale
            for (k = 0; k < params->box->n_tilings; k++)
            {
                for (l = 0; l < N_DIMS; l++)
                    temp[l] = positions[index * N_DIMS + l] + params->box->tilings[k * N_DIMS + l];
                unscale_coords(&positions[nc * N_DIMS * (k + 1) + index * N_DIMS],
                               temp, params->box);
            }

            // unscale the non-tiled folded scaled coordinates
            unscale_coords(temp,
                           &positions[index * N_DIMS], params->box);
            memcpy(&positions[index * N_DIMS], temp, N_DIMS * sizeof(SCALAR));
        }
    }
    return index_offset + group->n_gparticles * (group->size - 1);
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
        j = _fold_particles(params, group, positions, i, j);
        i += group->n_gparticles;
    }
}

void apply_constraints(run_params_t *params, SCALAR *positions, SCALAR *velocities)
{
    // Compute and add lagrange multipliers
    // implied by zeroth member of groups
    // typically identity, so irrelevant
    // except for Wyckoff positions
    unsigned int i, j, k, index = 0;
    double lambda, delta;
    SCALAR temp[N_DIMS];
    for (group_t *group = params->box->group; group != NULL; group = group->next)
    {
#pragma omp parallel for default(shared) private(i, j, k, temp, lambda, delta, index)
        for (i = 0; i < group->n_gparticles; i++)
        {
            scale_wrap_coords(&params->scaled_positions[index * N_DIMS], &positions[index * N_DIMS], params->box);
            action(group->members[0].g, temp,
                   &params->scaled_positions[index * N_DIMS], N_DIMS, 1.0);
            for (j = 0; j < N_DIMS; j++)
                temp[j] -= params->scaled_positions[index * N_DIMS + j];
            for (j = 0; j < N_DIMS; j++)
            {
                lambda = 0;
                for (k = 0; k < N_DIMS; k++)
                    lambda += group->ainv[j * N_DIMS + k] * temp[k];
#ifdef DEBUG
                printf("Constraint force for %d (dim %d) in group %s: %f. Ainv = %f. Delta = %f\n", index, j, group->name, lambda, group->ainv[j * N_DIMS],
                       temp[j]);
#endif
                // term of m / delta T cancels out
                // add constraint force to velocity
                velocities[index * N_DIMS + j] -= lambda / params->time_step * 2;
#ifdef DEBUG
                printf("position before constraints: %f %f\n", positions[index * N_DIMS + 0], positions[index * N_DIMS + 1]);
#endif
                // add to positions
                positions[index * N_DIMS + j] -= lambda;
#ifdef DEBUG
                printf("position after constraints: %f %f\n", positions[index * N_DIMS + 0], positions[index * N_DIMS + 1]);
#endif
            }
            scale_wrap_coords(&params->scaled_positions[index * N_DIMS], &positions[index * N_DIMS], params->box);
            index++;
        }
    }
}

void update_group(group_t *group, box_t *box)
{
    // update lagrangian matrix for group
    // index meanings
    // i - euclidean coordinate
    // j - scaled coordinate
    // k - constraint
    // matrix is (p_kj * b_ji - b_ki)^2
    // b = inverse box vectors transposed
    unsigned int i, j, k, l;
    double a;

    gsl_matrix *grad = gsl_matrix_alloc(N_DIMS, N_DIMS);
    gsl_matrix *mat = gsl_matrix_alloc(N_DIMS, N_DIMS);

    for (i = 0; i < N_DIMS; i++)
        for (j = 0; j < N_DIMS; j++)
            gsl_matrix_set(mat, i, j, group->members[0].g[i * (N_DIMS + 1) + j]);

    memcpy(grad->data, box->ib_vectors, sizeof(SCALAR) * N_DIMS * N_DIMS);
    // gsl_matrix_transpose(grad);

    // grad = p @ ib.T - ib.T
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mat, grad, -1.0, grad);
    // mat = grad @ grad
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, grad, grad, 0.0, mat);

// print A matrix
#ifdef DEBUG

    printf("A for group %s:\n", group->name);
    for (i = 0; i < N_DIMS; i++)
    {
        for (j = 0; j < N_DIMS; j++)
        {
            printf("%f ", gsl_matrix_get(mat, i, j));
        }
        printf("\n");
    }
#endif

    // matrix inverse never works here - maybe I'm doing something wrong
    // so use SVD

    gsl_matrix *V = gsl_matrix_alloc(N_DIMS, N_DIMS);
    gsl_vector *S = gsl_vector_alloc(N_DIMS);
    gsl_vector *work = gsl_vector_alloc(N_DIMS);
    gsl_matrix *Sp = gsl_matrix_calloc(N_DIMS, N_DIMS);
    // mat = U
    // V is untransposed
    gsl_linalg_SV_decomp(mat, V, S, work);

    // discard singular values and compute inverse S
    a = 0;
    for (i = 0; i < N_DIMS; i++)
    {
        a = fmax(a, gsl_vector_get(S, i));
    }
    // results do not look good with this implementation
    // a *= N_DIMS * 2 * FLT_MIN;
    a *= 1e-15 * N_DIMS;
    for (i = 0; i < N_DIMS; i++)
    {
        if (gsl_vector_get(S, i) > a)
        {
            gsl_matrix_set(Sp, i, i, 1.0 / gsl_vector_get(S, i));
        }
        else
        {
            gsl_matrix_set(Sp, i, i, 0.0);
        }
    }
    // compute inverse
    gsl_matrix *inv = gsl_matrix_alloc(N_DIMS, N_DIMS);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Sp, mat, 0.0, inv);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, inv, 0.0, mat);

// print A inv matrix
#ifdef DEBUG
    printf("A inv for group %s:\n", group->name);
    for (i = 0; i < N_DIMS; i++)
    {
        for (j = 0; j < N_DIMS; j++)
        {
            printf("%f ", gsl_matrix_get(mat, i, j));
        }
        printf("\n");
    }
#endif
    // now multiply by grad term
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mat, grad, 0.0, inv);

// print A inv matrix
#ifdef DEBUG

    printf("A inv w/grad for group %s:\n", group->name);
    for (i = 0; i < N_DIMS; i++)
    {
        for (j = 0; j < N_DIMS; j++)
        {
            printf("%f ", gsl_matrix_get(inv, i, j));
        }
        printf("\n");
    }
#endif
    // copy result
    memcpy(group->ainv, inv->data, sizeof(double) * N_DIMS * N_DIMS);

    gsl_matrix_free(mat);
    gsl_matrix_free(inv);
    gsl_matrix_free(V);
    gsl_matrix_free(Sp);
    gsl_matrix_free(grad);
    gsl_vector_free(S);
    gsl_vector_free(work);

    if (group->next)
        update_group(group->next, box);
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
