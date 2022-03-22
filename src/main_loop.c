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

  printf("Info: You are running version %s of symd\n", VERSION);

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
  double *forces = calloc(N_DIMS * params->n_particles, sizeof(double));
  char xyz_file_comment[100];
  double penergy = 0;
  double kenergy = 0;
  double insta_temperature;
  double therm_conserved = 0;
  SCALAR *output = (SCALAR *)malloc(sizeof(SCALAR) * N_DIMS * params->n_cell_particles);

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

    // output position
    if (do_exit || i % params->position_log_period == 0)
    {
      sprintf(xyz_file_comment, "Frame: %d", i);
      log_xyz(params->positions_file, positions, xyz_file_comment, N_DIMS,
              params->n_particles + params->n_ghost_particles, params->n_particles + params->n_ghost_particles, 0);
    }

    if (do_exit)
    {
      fprintf(stderr, "Exiting due to high temperature or energy  (T = %g)\n", insta_temperature);
      exit(1);
    }

    if (i % params->velocity_log_period == 0)
      log_array(params->velocities_file, velocities, N_DIMS, params->n_particles + params->n_ghost_particles, true);

    // integrate B1
    integrate_vel(params->time_step, velocities, forces, params->masses, params->n_particles);

    // integrate A1
    integrate_pos(params->time_step, positions, velocities, params->n_particles);

    // apply constraints
    if (params->box->group)
      apply_constraints(params, positions, velocities);

    // integrate O
    // thermostat
    if (params->thermostat_parameters)
      therm_conserved += params->thermostat_parameters->thermo_fxn(params->thermostat_parameters, params->temperature, params->time_step, positions, velocities, params->masses, N_DIMS, params->n_particles);

    // integrate A2
    integrate_pos(params->time_step, positions, velocities, params->n_particles);

    // apply constraints
    if (params->box->group)
      apply_constraints(params, positions, velocities);

    // apply group fold
    if (params->box->group)
      fold_particles(params, positions);

    // apply NPT step (before forces)
    if (params->box_update_period > 0 && i % params->box_update_period == 0)
      try_rescale(params, positions, &penergy, forces);

    // integrate B2
    // gather forces
    if (params->force_parameters)
      penergy = params->force_parameters->gather(params, positions, forces);
    integrate_vel(params->time_step, velocities, forces, params->masses, params->n_particles);

    // output forces
    if (i % params->force_log_period == 0)
      log_array(params->forces_file, forces, N_DIMS, params->n_particles + params->n_ghost_particles, true);

    // calculate important quantities
    kenergy = calculate_kenergy(velocities, params->masses, N_DIMS, params->n_particles);
    insta_temperature = kenergy * 2 / params->dof;

    if (i % params->print_period == 0)
    {
      printf("%12d %12g %12g %12g %12g %12g %12g %12g\n",
             i, i * params->time_step, insta_temperature, penergy, kenergy,
             penergy + kenergy, penergy + kenergy - therm_conserved, volume(params->box));

      // store "final" values at print period in case of bad exit
      for (j = 0; j < params->n_cell_particles; j++)
        scale_wrap_coords(&output[j * N_DIMS], &positions[j * N_DIMS], params->box);

      fseek(params->final_positions_file, 0, SEEK_SET);
      fseek(params->cell_file, 0, SEEK_SET);
      log_array(params->final_positions_file, output, N_DIMS,
                params->n_cell_particles, false);

      log_array(params->cell_file, params->box->b_vectors, N_DIMS,
                N_DIMS, false);
      log_array(params->cell_file, params->box->unorm_b_vectors, N_DIMS,
                N_DIMS, false);
    }
    if (params->n_particles > 1 && (insta_temperature != insta_temperature || insta_temperature > 1e25))
    {
      do_exit = 1;
    }
  }

  // only store if we had a clean exit
  if (!do_exit)
  {
    fseek(params->final_positions_file, 0, SEEK_SET);
    fseek(params->cell_file, 0, SEEK_SET);
    // want to store whole cell
    for (i = 0; i < params->n_cell_particles; i++)
      scale_wrap_coords(&output[i * N_DIMS], &positions[i * N_DIMS], params->box);
    log_array(params->final_positions_file, output, N_DIMS,
              params->n_cell_particles, false);

    log_array(params->cell_file, params->box->b_vectors, N_DIMS,
              N_DIMS, false);
    log_array(params->cell_file, params->box->unorm_b_vectors, N_DIMS,
              N_DIMS, false);
  }

  free(forces);
  free(output);
}
