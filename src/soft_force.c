#include "force.h"
#include <math.h>

#define PI 3.14159265

/* A soft cosine potential -- E = cos(pi * r).
 * F = -sin(pi * r) * x / r
 *
 */

double soft_gather_forces(run_params_t *params, double *positions, double *forces)
{

	unsigned int n_dims = params->n_dims, n_particles = params->n_particles;
	unsigned int n_ghost_particles = params->n_ghost_particles;
	double *box_size = params->box->box_size;
	unsigned int i, j, k;
	double penergy = 0;
	double r;
	double force_vector[n_dims];

	//zero forces
#pragma omp parallel for
	for (i = 0; i < n_particles; i++)
		for (k = 0; k < n_dims; k++)
			forces[i * n_dims + k] = 0;

#pragma parallel for default(shared) private(i, j, r, k, force_vector) \
	reduction(+                                                        \
			  : penergy)
	for (i = 0; i < n_particles; i++)
	{
		for (j = i + 1; j < n_particles + n_ghost_particles; j++)
		{
			r = 0;
			for (k = 0; k < n_dims; k++)
			{
				force_vector[k] = min_image_dist(positions[j * n_dims + k] - positions[i * n_dims + k], box_size[k]);
				//force_vector[k] = positions[j * n_dims + k] - positions[i * n_dims + k];
				r += force_vector[k] * force_vector[k];
			}
			r = sqrt(r);
			if (r < 0.00000001)
			{ //perturb
				force_vector[0] += 0.01;
				r = 0.01;
			}
			if (r < 1.5)
			{
#pragma omp critical
				for (k = 0; k < n_dims; k++)
				{
					forces[i * n_dims + k] += -sin(PI * r / 1.5) * PI * force_vector[k] / r;
					if (j < params->n_particles)
						forces[j * n_dims + k] += sin(PI * r / 1.5) * PI * force_vector[k] / r;
				}
#ifdef DEBUG
				printf("F(%d <-> %d, %g) = %g\n", i, j, r, -sin(PI * r / 1.5));
#endif //DEBUG
				penergy += cos(PI * r / 1.5);
			}
			else
				penergy += 0;
		}
	}

	return (penergy);
}

void s_free_forces(force_t *parameters)
{
	return;
}

force_t *build_soft()
{
	force_t *soft = (force_t *)malloc(sizeof(force_t));
	soft->gather = &soft_gather_forces;
	soft->free = &s_free_forces;
	soft->parameters = NULL;
}