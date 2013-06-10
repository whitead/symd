#include "main_loop.h"
#include <stdio.h>

int main(int argc, char* argv[]) {

  char* pfile;
  //process arguemnts to get parameters
  pfile = NULL;

  if(argc == 1)
    pfile = argv[1];

  Run_Params* p = read_parameters(pfile);

  //start main loop
  main_loop(p);

  free_run_params(p);

  return(0);

}

int main_loop(Run_Params* params){

  unsigned int i;
  double* positions = params->initial_positions;
  double* velocities = params->initial_velocities;
  double* forces = malloc(sizeof(double) * params->n_dims * params->n_particles);
  char xyz_file_comment[100];
  double penergy = 0;
  double kenergy = 0;
  double insta_temperature;
  double therm_conserved = 0;

  gather_forces(params->force_parameters, positions, forces, params->masses, params->box_size, params->n_dims, params->n_particles);


  printf("%12s %12s %12s %12s %12s %12s %12s\n", "Step", "Time", "T", "PE", "KE", "E", "Htherm");
  //start at 0, so that we don't log on the first loop
  for(i = 1; i <= params->steps; i++) {

    //integrate 1
    integrate_1(params->time_step, positions, velocities, forces, params->masses,params->box_size,  params->n_dims, params->n_particles);

    //remove COM if necessary
    if(i % params->com_remove_period == 0)
      remove_com(velocities, params->masses, params->n_dims, params->n_particles);

    //gather forces
    penergy =  gather_forces(params->force_parameters, positions, forces, params->masses, params->box_size,  params->n_dims, params->n_particles);

    //integrate 2
    integrate_2(params->time_step, positions, velocities, forces, params->masses, params->box_size, params->n_dims, params->n_particles);

    //thermostat
#ifdef THERMOSTAT
    therm_conserved += thermostat(params->temperature, params->time_step, params->thermostat_parameters, positions, velocities, params->masses, params->n_dims, params->n_particles);
#endif

    //calculate important quantities
    kenergy = calculate_kenergy(velocities, params->masses, params->n_dims, params->n_particles);
    insta_temperature = kenergy * 2 / (params->n_particles * params->n_dims - params->n_dims);
    
    
    //print
    if(i % params->position_log_period == 0) {      
      sprintf(xyz_file_comment, "Frame: %d", i);
	log_xyz(params->positions_file, positions, xyz_file_comment, params->n_dims, params->n_particles);
    }

    if(i % params->velocity_log_period == 0)
      log_array(params->velocities_file, velocities, params->n_dims, params->n_particles, true);

    if(i % params->force_log_period == 0)
      log_array(params->forces_file, forces, params->n_dims, params->n_particles, true);

    if(i % params->print_period == 0) {
      printf("%12d %12g %12g %12g %12g %12g %12g\n", i, i * params->time_step, insta_temperature, penergy, kenergy, penergy + kenergy, penergy + kenergy - therm_conserved);
    }
  }

  if(params->positions_file != NULL)
    fclose(params->positions_file);
  if(params->velocities_file != NULL)
    fclose(params->velocities_file);
  if(params->forces_file != NULL)
    fclose(params->forces_file);

  free(forces);
  
  return MD_SUCCESS;
}


