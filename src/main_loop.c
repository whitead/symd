#include "main_loop.h"
#include <stdio.h>
#include <stdlib.h>
#include "force.h"
#include "integrate.h"
#include "thermostat.h"
#include "group.h"
#include "box.h"

int main(int argc, char *argv[])
{

  char *pfile;
  // process arguments to get parameters
  pfile = NULL;

  if (argc == 2)
    pfile = argv[1];

  printf("You are running version %s of symd\n", VERSION);

  run_params_t *p = read_parameters(pfile);

  // start main loop
  main_loop(p);

  free_run_params(p);

  return (0);
}

void main_loop(run_params_t *params)
{

  unsigned int i, j, k;
  int do_exit = 0;
  double *positions = params->initial_positions;
  double *velocities = params->initial_velocities;
  // make sure initial forces are zeroed
  double *forces = calloc(params->n_dims * params->n_particles, sizeof(double));
  char xyz_file_comment[100];
  double penergy = 0;
  double kenergy = 0;
  double insta_temperature;
  double therm_conserved = 0;

  // apply group if necessary
  if (params->box->group)
    fold_particles(params, positions);

  // gather forces -> must be here so we don't do integrate 1 without forces
  if (params->force_parameters)
    penergy = params->force_parameters->gather(params, positions, forces);

  printf("%12s %12s %12s %12s %12s %12s %12s %12s\n",
         "Step", "Time", "T", "PE", "KE", "E", "Htherm",
         "V");
  // start at 0, so that we don't log on the first loop
  for (i = 0; i < params->steps; i++)
  {

    // output
    if (do_exit || i % params->position_log_period == 0)
    {
      sprintf(xyz_file_comment, "Frame: %d", i);
      for (j = 0; j < params->box->group->size; j++)
      {
        // asymmetric
        if (j == 0)
          log_xyz(params->positions_file, positions, xyz_file_comment, elements, params->n_dims,
                  params->n_particles, params->n_particles + params->n_ghost_particles, 0);
        else // non-tiled
          log_xyz(params->positions_file, &positions[params->n_dims * params->n_particles * j], NULL, elements, params->n_dims,
                  params->n_particles, params->n_particles + params->n_ghost_particles, 1);
        for (k = 0; k < params->box->n_tilings; k++)
          log_xyz(params->positions_file,
                  &positions[params->n_dims * (params->n_particles * ((k + 1) * params->box->group->size + j))],
                  NULL, elements,
                  params->n_dims, params->n_particles,
                  params->n_particles + params->n_ghost_particles, 1);
      }
    }

    if (do_exit)
    {
      fprintf(stderr, "Exiting due to high temperature or energy\n");
      exit(1);
    }

    if (i % params->velocity_log_period == 0)
      log_array(params->velocities_file, velocities, params->n_dims, params->n_particles + params->n_ghost_particles, true);

    // integrate 1
    integrate_1(params->time_step, positions, velocities, forces, params->masses, params->box->box_size, params->n_dims, params->n_particles);

    // apply group if necessary
    if (params->box->group)
      fold_particles(params, positions);

    if (i % params->com_remove_period == 0)
      remove_com(velocities, params->masses, params->n_dims, params->n_particles);

    // gather forces
    if (params->force_parameters)
      penergy = params->force_parameters->gather(params, positions, forces);

    // try update box
    // if (params->box_update_period && i % params->box_update_period == 0)
    // {
    //   if (!try_rescale(params, positions, &penergy, forces))
    //   {
    //     // reset forces
    //     penergy = params->force_parameters->gather(params, positions, forces);
    //   }
    // }

    // output forces
    if (i % params->force_log_period == 0)
      log_array(params->forces_file, forces, params->n_dims, params->n_particles + params->n_ghost_particles, true);

    // integrate 2
    integrate_2(params->time_step, positions, velocities, forces, params->masses, params->box->box_size, params->n_dims, params->n_particles);

    // thermostat
    if (params->thermostat_parameters)
      therm_conserved += params->thermostat_parameters->thermo_fxn(params->thermostat_parameters, params->temperature, params->time_step, positions, velocities, params->masses, params->n_dims, params->n_particles);

    // calculate important quantities
    kenergy = calculate_kenergy(velocities, params->masses, params->n_dims, params->n_particles);
    insta_temperature = kenergy * 2 / (params->n_particles * params->n_dims - params->n_dims);

    if (i % params->print_period == 0)
    {
      printf("%12d %12g %12g %12g %12g %12g %12g %12g\n",
             i, i * params->time_step, insta_temperature, penergy, kenergy,
             penergy + kenergy, penergy + kenergy - therm_conserved, volume(params->box->box_size, params->n_dims));
    }
    if (insta_temperature > 1000)
    {
      do_exit = 1;
    }
  }
  log_array(params->final_positions_file, positions, params->n_dims,
            params->n_particles + params->n_ghost_particles, false);

  free(forces);
}
