
#include <stdio.h>

#ifndef VERSION
#define VERSION unknown
#endif

#ifndef PARAMS_H_
#define PARAMS_H_

typedef struct force_t force_t;
typedef struct thermostat_t thermostat_t;
typedef struct group_t group_t;

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
    thermostat_t *thermostat_parameters;
    force_t *force_parameters;
    group_t *group;
    double *initial_positions;
    double *initial_velocities;
    double *masses;
    double *box_size;
    unsigned int n_dims;
    unsigned int n_particles;
    unsigned int n_ghost_particles;
    FILE *positions_file;
    FILE *final_positions_file;
    FILE *velocities_file;
    FILE *forces_file;
    unsigned int position_log_period;
    unsigned int velocity_log_period;
    unsigned int force_log_period;
    unsigned int print_period;

} run_params_t;

#endif