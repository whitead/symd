#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include "thermostat.h"
#include "cJSON.h"
#include "params.h"

#ifndef UTIL_H_
#define UTIL_H_

/*
 * Generate a set of initial velocities from the appropiate chi-n distribution
 *
 */
double *generate_velocities(double temperature, unsigned int seed, double *masses, unsigned int n_dims, unsigned int n_particles);

/*
 * Calculate the unitless kinetic energy.
 *
 */
double calculate_kenergy(double *velocities, double *masses, unsigned int n_dims, unsigned int n_particles);

/*
 * Write to file, which should be open, the given array in xyz format
 */
void log_xyz(FILE *file, double *array, char *frame_string, unsigned n_dims, unsigned n_particles);

/*
 * Write to file, which should be open, the given array. If sum is true, then the cols will
 * be squared, summed, and the square root taken. For example, it will produce the speed for a
 * velocity array.
 */
void log_array(FILE *file, double *array, unsigned n_cols, unsigned n_rows, bool sum);

/*
 * Load a matrix. Skip is the number of lines to skip before beginning.
 *
 */
double *load_matrix(char *filename, unsigned int nrow, unsigned int ncol, unsigned int skip);

/*
 * Read in the given parameters.
 */
run_params_t *read_parameters(char *params_file);

/*
* Load group from JSON file
*/
group_t *load_group(char *filename, unsigned int n_dims);

unsigned int process_uint(char ***pstrings, char *key, bool *success);

double process_double(char ***pstrings, char *key, bool *success);

char *process_string(char ***pstrings, char *key, bool *success);

/* Remove center of mass translational motion. Returns the magnitude of the COM
 *
 */
double remove_com(double *velocities, double *masses, unsigned int n_dims, unsigned int n_particles);

//free all the memory bits in the run_params_t struct
void free_run_params(run_params_t *params);

#endif
