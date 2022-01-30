#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include "thermostat.h"
#include "force.h"
#include "cJSON.h"

#ifndef UTIL_H_
#define UTIL_H_

#define MD_SUCCESS 0

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
  double *initial_positions;
  double *initial_velocities;
  double *masses;
  double *box_size;
  unsigned int n_dims;
  unsigned int n_particles;
  FILE *positions_file;
  FILE *velocities_file;
  FILE *forces_file;
  unsigned int position_log_period;
  unsigned int velocity_log_period;
  unsigned int force_log_period;
  unsigned int print_period;

} Run_Params;

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
Run_Params *read_parameters(char *params_file);

unsigned int process_uint(char ***pstrings, char *key, bool *success);

double process_double(char ***pstrings, char *key, bool *success);

char *process_string(char ***pstrings, char *key, bool *success);

/* Remove center of mass translational motion. Returns the magnitude of the COM
 *
 */
double remove_com(double *velocities, double *masses, unsigned int n_dims, unsigned int n_particles);

//free all the memory bits in the Run_Params struct
void free_run_params(Run_Params *params);

#endif
