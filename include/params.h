
#include <stdio.h>
#include <gsl/gsl_rng.h>

#ifndef VERSION
#define VERSION unknown
#endif

#ifndef PARAMS_H_
#define PARAMS_H_

typedef struct force_t force_t;
typedef struct thermostat_t thermostat_t;
typedef struct group_t group_t;
typedef struct box_t box_t;

typedef double SCALAR;

/*
 * Parameters needed to begin an MD run
 *
 */
typedef struct
{

    unsigned int steps;
    unsigned int com_remove_period;
    double time_step;
    double temperature;
    double pressure;
    thermostat_t *thermostat_parameters;
    force_t *force_parameters;
    box_t *box;
    double *initial_positions;
    double *initial_velocities;
    double *scaled_positions;
    double *masses;
    unsigned int n_particles;
    unsigned int n_ghost_particles;
    FILE *positions_file;
    FILE *final_positions_file;
    FILE *velocities_file;
    FILE *forces_file;
    unsigned int position_log_period;
    unsigned int velocity_log_period;
    unsigned int force_log_period;
    unsigned int box_update_period;
    unsigned int print_period;
    gsl_rng *rng;

} run_params_t;

#endif