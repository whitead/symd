#include "util.h"

#define ARG_BUFFER 150
#define STRING_BUFFER 50
#define KEY 0
#define VALUE 1

double* generate_velocities(double temperature, unsigned int seed, double* masses, unsigned int n_dims, unsigned int n_particles) {

  double* velocities = (double*) malloc(sizeof(double) * n_dims * n_particles);

  unsigned int i, j;
  const gsl_rng * rng;

  gsl_rng_env_setup();
  rng = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(rng, seed);
  
  for(i = 0; i < n_particles; i++) {
    for(j = 0; j < n_dims; j++) {
      velocities[i * n_dims + j] = gsl_ran_gaussian_ziggurat(rng, sqrt(temperature / masses[i]));
    }
  }

  double ke = calculate_kenergy(velocities, masses, n_dims, n_particles);

  printf("Generated velocity distribution with %g temperature\n", ke * 2 / (n_dims * n_particles));

  return(velocities);

}


Run_Params* read_parameters(FILE* params_file, const Run_Params* default_params) {


  //set up a buffer for all the arguments. The first dimension is file row, the second
  // is 2, one for key and one for value. The final dimension is the string itself.
  char*** pstrings = (char***) malloc(sizeof(char**) * ARG_BUFFER);
  unsigned int i,j;

  for(i = 0; i < ARG_BUFFER; i++) {
    pstrings[i] = (char**) malloc(sizeof(char *) * 2);
    for(j = 0; j < 2; j++) {
      pstrings[i][j] = (char*) malloc(sizeof(char) * STRING_BUFFER);
    }
  }

  unsigned int n = 2;

  for(i = 0; n == 2 && i < ARG_BUFFER; i++) {
    n = fscanf(params_file,"%s %s\n", pstrings[i][KEY], pstrings[i][VALUE]);
    printf("Read [%s] and [%s]\n", pstrings[i][KEY], pstrings[i][VALUE]);
  }

  Run_Params* params = (Run_Params*) malloc(sizeof(Run_Params));

  //Set to default if non-null
  if(default_params != NULL) {
    memcpy(params, default_params, sizeof(Run_Params));
  }

  //Used to keep track of success in finding key/values
  bool success_cur;
  
  //store in temps and only store if successfully found value in param file
  //This allows the default value to not be overwritten if things were unsuccessful
  unsigned int temp_uint;
  double temp_double;

  //Get unsigned integers

  temp_uint = process_uint(pstrings, "steps", &success_cur);
  if(success_cur)
    params->steps = temp_uint;
  else if(default_params == NULL) {
    fprintf(stderr, "Could not read steps\n");
    exit(1);
  }

  temp_uint = process_uint(pstrings, "n_dims", &success_cur);
  if(success_cur)
    params->n_dims = temp_uint;
  else if(default_params == NULL) {
    fprintf(stderr, "Could not read n_dims\n");
    exit(1);
  }

  temp_uint = process_uint(pstrings, "n_particles", &success_cur);
  if(success_cur)
    params->n_particles = temp_uint;
  else if(default_params == NULL) {
    fprintf(stderr, "Could not read n_particles\n");
    exit(1);
  }

  temp_uint = process_uint(pstrings, "print_period", &success_cur);
  if(success_cur)
    params->print_period = temp_uint;
  else if(default_params == NULL) {
    fprintf(stderr, "Could not read print_period\n");
    exit(1);
  }

  temp_uint = process_uint(pstrings, "position_log_period", &success_cur);
  if(success_cur)
    params->position_log_period = temp_uint;
  else if(default_params == NULL) {
    params->position_log_period = params->print_period;
  }


  temp_uint = process_uint(pstrings, "velocity_log_period", &success_cur);
  if(success_cur)
    params->velocity_log_period = temp_uint;
  else if(default_params == NULL) {
    params->velocity_log_period = params->print_period;
  }

  temp_uint = process_uint(pstrings, "force_log_period", &success_cur);
  if(success_cur)
    params->force_log_period = temp_uint;
  else if(default_params == NULL) {
    params->force_log_period = params->print_period;
  }


  temp_double = process_double(pstrings, "time_step", &success_cur);
  if(success_cur)
    params->time_step = temp_double;
  else if(default_params == NULL) {
    fprintf(stderr, "Could not read time_step\n");
    exit(1);
  }

  temp_double = process_double(pstrings, "temperature", &success_cur);
  if(success_cur)
    params->temperature = temp_double;
  else if(default_params == NULL) {
    fprintf(stderr, "Could not read temperature\n");
    exit(1);
  }


  //thermostat parameters
  #ifdef ANDERSON
  unsigned int seed = process_uint(pstrings, "anderson_seed", &success_cur);
  if(!success_cur) {
    seed = 1523;
    fprintf(stderr, "Warning: Using default seed for anderson thermostat\n");
  }
  double collision_freq = process_double(pstrings, "anderson_nu", &success_cur);
  if(!success_cur) {
    collision_freq = 10;
  }
  
  params->thermostat_parameters = build_anderson(seed, collision_freq);
  #endif
  

  //System parameters
  #ifdef HARMONIC
  double k = process_double(pstrings, "harmonic_constant", &success_cur);
  if(!success_cur) {
    k = 1;
  }
  params->force_parameters = build_harmonic(k);
  #endif

  #ifdef LJ
  double epsilon = process_double(pstrings, "lj_epsilon", &success_cur);  
  double sigma = process_double(pstrings, "lj_sigma", &success_cur);
  if(!success_cur) {
    epsilon = 1;
    sigma = 1;
    fprintf(stderr, "Could not find LJ parameters. Assuming sigma = epsilon = 1\n");
  }
  params->force_parameters = build_lj(epsilon, sigma);
  #endif

  
  //load starting configuration

  char* positions_file = process_string(pstrings, "start_positions", &success_cur);
  if(!success_cur) {
    fprintf(stderr, "Could not find start_positions file name\n");
    exit(1);
  }
  
  params->initial_positions = load_matrix(positions_file, params->n_particles, params->n_dims, 0);  

  //load masses
  
  char* masses_file = process_string(pstrings, "masses_file", &success_cur);
  if(!success_cur) {
    fprintf(stderr, "Could not find masses_file name\n");
    exit(1);
  }

  params->masses = load_matrix(masses_file, params->n_particles, 1, 0);  
  
  //load, or create, velocities

  char* velocities_file = process_string(pstrings, "start_velocities", &success_cur);
  if(!success_cur) {
    temp_uint = process_uint(pstrings, "velocity_seed", &success_cur);
    if(!success_cur) {
      fprintf(stderr, "Warning: Using default velocity seed\n");
      temp_uint = 52989;
    }
    params->initial_velocities = 
      generate_velocities(params->temperature, temp_uint, 
			  params->masses, params->n_dims,
			  params->n_particles);
  } else{ 
    params->initial_velocities = 
      load_matrix(velocities_file, params->n_particles, params->n_dims, 0);
  }
  
  //read in output file names. If there are no values present, the process_string returns
  //NULL which is assumed to be no output desired.

  char* temp_string;



  temp_string = process_string(pstrings, "positions_log_file", &success_cur);
  if(temp_string != NULL) {
    params->positions_file = fopen(temp_string, "w");
  }
  else {
    params->positions_file = NULL;
  }

  temp_string = process_string(pstrings, "velocities_log_file", &success_cur);
  if(temp_string != NULL) {
    params->velocities_file = fopen(temp_string, "w");
  }
  else {
    params->velocities_file = NULL;
  }

  temp_string = process_string(pstrings, "forces_log_file", &success_cur);
  if(temp_string != NULL) {
    params->forces_file = fopen(temp_string, "w");
  }
  else {
    params->forces_file = NULL;
  }    


  for(i = 0; i < ARG_BUFFER; i++) {
    for(j = 0; j < 2; j++) {
      free(pstrings[i][j]);
      pstrings[i][j] = NULL;
    }
    free(pstrings[i]);
    pstrings[i] = NULL;
  }
  free(pstrings);
  pstrings = NULL;

  
  return params;
  
}

char* process_string(char*** pstrings, char* key, bool* success) {

  unsigned int i;
  for(i = 0; i < ARG_BUFFER; i++) {
    
    if(strcmp(pstrings[i][KEY], key) == 0) {
      char*  result = pstrings[i][VALUE];
      *success = true;
      printf("Read in [%s] for [%s]\n", result, key);
      return result;
    }
  }

  *success = false;

  return NULL;
}

unsigned int process_uint(char*** pstrings, char* key, bool* success) {

  unsigned int i;
  for(i = 0; i < ARG_BUFFER; i++) {
    
    if(strcmp(pstrings[i][KEY], key) == 0) {
      unsigned int result;
      sscanf(pstrings[i][VALUE], "%ud", &result);
      *success = true;
      printf("Read in [%d] for [%s]\n", result, key);
      return result;
    }
  }

  *success = false;

  return 0;
}

double process_double(char*** pstrings, char* key, bool* success) {

  unsigned int i;
  for(i = 0; i < ARG_BUFFER; i++) {
    
    if(strcmp(pstrings[i][KEY], key) == 0) {
      double result;
      sscanf(pstrings[i][VALUE], "%lf", &result);
      *success = true;
      printf("Read in [%6f] for [%s]\n", result, key);
      return result;
    }
  }

  *success = false;

  return 0;
}


double calculate_kenergy(double* velocities, double* masses, unsigned int n_dims, unsigned int n_particles) {

  unsigned int i, j;
  double kenergy = 0;
  double etemp;
  for(i = 0; i < n_particles; i++) {
    etemp = 0;
    for(j = 0; j < n_dims; j++) {
      etemp += velocities[i * n_dims + j] * velocities[i * n_dims + j];
    }
    kenergy += etemp * masses[i] / 2;
  }

  return(kenergy);

}

void log_xyz(FILE* file, double* array, char* frame_string, unsigned n_dims, unsigned n_particles) {
  
  if(file == NULL)
    return;

  unsigned int i,j;

  fprintf(file, "%d\n%s\n", n_particles, frame_string);

  for(i = 0; i < n_particles; i++) {
    fprintf(file, "Ar ");
    for(j = 0; j < n_dims; j++) {
      fprintf(file, "%12g ", array[i * n_dims + j]); 
    }
      fprintf(file, "\n");
  }

  fflush(file);
  
}

void log_array(FILE* file, double* array, unsigned n_cols, unsigned n_rows, bool do_sum) {
  
  if(file == NULL)
    return;

  unsigned int i,j;
  double sum;

  for(i = 0; i < n_rows; i++) {
    if(do_sum)
      sum = 0;
    for(j = 0; j < n_cols; j++) {
      if(do_sum)
	sum += array[i * n_cols + j] * array[i * n_cols + j];
      fprintf(file, "%12g ", array[i * n_cols + j]); 
    }
    if(do_sum)
      fprintf(file, "%12g\n", sqrt(sum));
    else
      fprintf(file, "\n");
  }

  fflush(file);
  
}

double* load_matrix(char* filename, unsigned int nrow, unsigned int ncol, unsigned int skip) {

  FILE *mfile;
  mfile = fopen(filename, "r");
  if(mfile != NULL) {

    int c;
    unsigned int i, j;

    if(skip > 0) {
      for(i = 0; i < skip; i++) {
	while((c = getc(mfile)) != '\n' && c != EOF);
      }
    }

    double *matrix =  (double*) malloc(sizeof(double) * nrow * ncol);
    if(matrix == NULL) {
      perror("Out of memory\n");
      exit(1);
    }

    for(i = 0; i < nrow; i++) {
      for(j = 0; j < ncol; j++) {
	if(fscanf(mfile, "%lg", &(matrix[i*ncol + j])) == 0) {
	  fprintf(stderr, "Incorrect number of rows or columns"
		  "at i = %d, and j=%d, nrow=%d, ncol=%d\n", i, j, nrow, ncol);
	  exit(1);
	}
      }  
    }

    fclose(mfile);
    return matrix;
    
  }else {
    perror("Could not open file\n");
  }

  return NULL;

}


void free_run_params(Run_Params* params) {

  free(params->thermostat_parameters);
  free(params->force_parameters);
  free(params->initial_positions);
  free(params->initial_velocities);
  free(params->masses);

  free(params);

}
