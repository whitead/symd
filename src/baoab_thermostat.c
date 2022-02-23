#include "thermostat.h"
#include <gsl/gsl_randist.h>
#include <math.h>

double baoab_thermostat(thermostat_t *params, double temperature, double time_step, double *positions, double *velocities, double *masses, unsigned int n_dims, unsigned int n_particles)
{

    unsigned int i, j;
    double c1 = exp(-params->param * time_step);
    double c2 = sqrt(1 - c1 * c1);
    for (i = 0; i < n_particles; i++)
    {
        for (j = 0; j < n_dims; j++)
        {
            velocities[i * n_dims + j] *= c1;
            velocities[i * n_dims + j] += c2 * gsl_ran_gaussian_ziggurat(params->rng, sqrt(temperature / masses[i]));
        }
    }

    return 0;
}

void baoab_free_thermostat(thermostat_t *params)
{
    free(params);
}

thermostat_t *build_baoab(double gamma, gsl_rng *rng)
{
    thermostat_t *params = (thermostat_t *)malloc(sizeof(thermostat_t));
    params->param = gamma;
    params->rng = rng;
    params->thermo_fxn = &baoab_thermostat;
    params->free = &baoab_free_thermostat;
    return params;
}
