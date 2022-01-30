#include "main_loop.h"
#include <stdio.h>
#include <stdlib.h>
#include "force.h"
#include "integrate.h"
#include "thermostat.h"
#include "group.h"

int main(int argc, char *argv[])
{

  char *pfile;
  //process arguemnts to get parameters
  pfile = NULL;

  if (argc == 2)
    pfile = argv[1];

  run_params_t *p = read_parameters(pfile);

  //start main loop
  main_loop(p);

  free_run_params(p);

  return (0);
}

int main_loop(run_params_t *params)
{

  unsigned int i;
  double *positions = params->initial_positions;
  double *velocities = params->initial_velocities;
  double *forces = malloc(sizeof(double) * params->n_dims * params->n_particles);
  char xyz_file_comment[100];
  double penergy = 0;
  double kenergy = 0;
  double insta_temperature;
  double therm_conserved = 0;

  // apply group if necessary
  if (params->group)
  {
    fold_particles(params, positions, true);
  }

  printf("%12s %12s %12s %12s %12s %12s %12s\n", "Step", "Time", "T", "PE", "KE", "E", "Htherm");
  //start at 0, so that we don't log on the first loop
  for (i = 1; i <= params->steps; i++)
  {

    //integrate 1
    integrate_1(params->time_step, positions, velocities, forces, params->masses, params->box_size, params->n_dims, params->n_particles);

    // apply group if necessary
    if (params->group)
    {
      fold_particles(params, positions, false);
      //fold_particles(params, velocities, false);
    }
    if (i % params->com_remove_period == 0)
      remove_com(velocities, params->masses, params->n_dims, params->n_particles);

    //gather forces
    penergy = params->force_parameters->gather(params, positions, forces);

    //integrate 2
    integrate_2(params->time_step, positions, velocities, forces, params->masses, params->box_size, params->n_dims, params->n_particles);

    //thermostat
    if (params->thermostat_parameters)
      therm_conserved += params->thermostat_parameters->thermo_fxn(params->thermostat_parameters, params->temperature, params->time_step, positions, velocities, params->masses, params->n_dims, params->n_particles);

    //calculate important quantities
    kenergy = calculate_kenergy(velocities, params->masses, params->n_dims, params->n_particles);
    insta_temperature = kenergy * 2 / (params->n_particles * params->n_dims - params->n_dims);

    //print
    if (i % params->position_log_period == 0)
    {
      sprintf(xyz_file_comment, "Frame: %d", i);
      log_xyz(params->positions_file, positions, xyz_file_comment, params->n_dims, params->n_particles + params->n_ghost_particles);
    }

    if (i % params->velocity_log_period == 0)
      log_array(params->velocities_file, velocities, params->n_dims, params->n_particles + params->n_ghost_particles, true);

    if (i % params->force_log_period == 0)
      log_array(params->forces_file, forces, params->n_dims, params->n_particles + params->n_ghost_particles, true);

    if (i % params->print_period == 0)
    {
      printf("%12d %12g %12g %12g %12g %12g %12g\n", i, i * params->time_step, insta_temperature, penergy, kenergy, penergy + kenergy, penergy + kenergy - therm_conserved);
    }
  }

  free(forces);

  return MD_SUCCESS;
}
