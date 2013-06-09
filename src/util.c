#include "util.h"


#define ARG_BUFFER 150
#define STRING_BUFFER 50
#define KEY 0
#define VALUE 1

static const char*  
default_json       = " { \"com_remove_period\" : 1000, \"skin\" : 0, \"thermostat_seed\" : 1523, \"anderson_nu\" : 10.0, \"harmonic_constant\" : 1.0, \"lj_epsilon\" : 1.0, \"lj_sigma\" : 1.0,  \"velocity_seed\" : 543214, \"position_log_period\" : 0, \"velocity_log_period\" : 0, \"force_log_period\" : 0} ";


double* generate_velocities(double temperature, unsigned int seed, double* masses, unsigned int n_dims, unsigned int n_particles) {

  double* velocities = (double*) malloc(sizeof(double) * n_dims * n_particles);

  unsigned int i, j;
  gsl_rng * rng;

  gsl_rng_env_setup();
  rng = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(rng, seed);
  
  for(i = 0; i < n_particles; i++) {
    for(j = 0; j < n_dims; j++) {
      if(temperature != 0)
	velocities[i * n_dims + j] = gsl_ran_gaussian_ziggurat(rng, sqrt(temperature / masses[i]));
      else
	velocities[i * n_dims + j] = 0;
    }
  }


#ifdef DEBUG
  double ke = calculate_kenergy(velocities, masses, n_dims, n_particles);
  printf("Generated velocity distribution with %g temperature\n", ke * 2 / (n_dims * n_particles));
#endif

  gsl_rng_free(rng);

  return(velocities);

}


void load_json(char* filename, char** data) {

  FILE* f;

  if(!filename) { //stdin    
    f = stdin;    
  } else {//file
    //get file length
    f = fopen(filename,"rb");
  }

  fseek(f,0,SEEK_END);
  long len = ftell(f);
  //reset file
  fseek(f,0,SEEK_SET);
  //load data
  *data = (char*) malloc(len+1);
  fread(*data,1,len,f);

  if(filename)
    fclose(f);

}

cJSON* retrieve_item(cJSON* root, cJSON* default_root, const char* item_name) {
  cJSON* item = cJSON_GetObjectItem(root, item_name);
  if(!item) { //check if it's in the default parameter list
    item = cJSON_GetObjectItem(default_root, item_name);
    if(item) {
      fprintf(stderr, "Warning: assuming default value for %s = ", item_name);
      switch(item->type) {
      case cJSON_False:
	fprintf(stderr, "false\n");
	break;
      case cJSON_True:
	fprintf(stderr, "true\n");
	break;
      case cJSON_Number:
	fprintf(stderr, "%g\n", item->valuedouble);
	break;
      case cJSON_String:
	fprintf(stderr, "%s\n", item->valuestring);
	break;
      }
	
    }
    else {
      fprintf(stderr, "Error: could not read %s and no default value\n", item_name);
      exit(1);
    }
  }

  return item;
}

Run_Params* read_parameters(char* file_name) {

  char* data;
  load_json(file_name, &data);
  cJSON* root = cJSON_Parse(data);
  cJSON* default_root = cJSON_Parse(default_json);
  cJSON* item;
  
  if (!root) {
    printf("Error before: [%s]\n",cJSON_GetErrorPtr());
    return NULL;
  }

  Run_Params* params = (Run_Params*) malloc(sizeof(Run_Params));


  //get single value parameters
  params->steps = (unsigned int) retrieve_item(root, default_root, "steps")->valueint;
  params->com_remove_period = (unsigned int) retrieve_item(root, default_root, "com_remove_period")->valueint;
  params->time_step =  retrieve_item(root, default_root, "time_step")->valuedouble;

  params->temperature = retrieve_item(root, default_root, "temperature")->valuedouble;

  params->n_dims = (unsigned int) retrieve_item(root, default_root, "n_dims")->valueint;
  params->n_particles = (unsigned int) retrieve_item(root, default_root, "n_particles")->valueint;

  params->print_period =  (unsigned int) retrieve_item(root, default_root, "print_period")->valueint;

  params->position_log_period =  (unsigned int) retrieve_item(root, default_root, "position_log_period")->valueint;
  params->position_log_period = (params->position_log_period == 0 ? params->print_period : params->position_log_period);

  params->velocity_log_period =  (unsigned int) retrieve_item(root, default_root, "velocity_log_period")->valueint;
  params->velocity_log_period = (params->velocity_log_period == 0 ? params->print_period : params->velocity_log_period);

  params->force_log_period =  (unsigned int) retrieve_item(root, default_root, "force_log_period")->valueint;
  params->force_log_period = (params->force_log_period == 0 ? params->print_period : params->force_log_period);

  
  //box size
#ifndef NO_PBC

  unsigned int i;
  params->box_size = (double*) malloc(sizeof(double) * params->n_dims);
  
  item = retrieve_item(root, default_root, "box_size");
  i = 0;
  for(item = item->child; item != NULL; item = item->next) {
    params->box_size[i] = item->valuedouble;
    i++;
  }
  if(i != params->n_dims) {
    fprintf(stderr, "Error: Number of box dimensions not equal to simulation dimension\n");
    exit(1);
  }
  
#endif

  //neighbor list parameters
#ifdef NLIST
  
  double rcut = retrieve_item(root, default_root, "rcut")->valuedouble;
  double skin = retrieve_item(root, default_root, "skin")->valuedouble;
  if(skin == 0) {
    skin = 0.2 * rcut;
    fprintf(stderr, "Warning: Assuming skin = %g\n", skin);
  }

  Nlist_Parameters* nlist = build_nlist_params(params->n_dims, params->n_particles,
					       params->box_size, skin, rcut);
#endif//NLIST

  //thermostats
#ifdef ANDERSON

  unsigned int seed = (unsigned int) retrieve_item(root, default_root, "thermostat_seed")->valueint;
  double collision_freq = retrieve_item(root, default_root, "anderson_nu")->valuedouble;
  params->thermostat_parameters = build_anderson(seed, collision_freq);

#endif//ANDERSON

#ifdef BUSSI

  unsigned int seed = (unsigned int) retrieve_item(root, default_root, "thermostat_seed")->valueint;
  double taut = retrieve_item(root, default_root, "bussi_taut")->valuedouble;
  params->thermostat_parameters = build_bussi(seed, taut);

#endif//BUSSI
  
#ifdef HARMONIC
  double k = retrieve_item(root, default_root, "harmonic_constant")->valuedouble;
  params->force_parameters = build_harmonic(k);
#endif//HARMONIC

#ifdef LJ
  double epsilon = retrieve_item(root, default_root, "lj_epsilon")->valuedouble;
  double sigma = retrieve_item(root, default_root, "lj_sigma")->valuedouble;
  params->force_parameters = build_lj(epsilon, sigma, nlist);
#endif//LJ

  //load input files
  char* positions_file = retrieve_item(root, default_root, "start_positions")->valuestring; 
  params->initial_positions = load_matrix(positions_file, params->n_particles, params->n_dims, 0);
  
  char* masses_file = retrieve_item(root, default_root, "masses_file")->valuestring;
  params->masses = load_matrix(masses_file, params->n_particles, 1, 0);

  item = cJSON_GetObjectItem(root, "start_velocities");
  if(!item) {
    unsigned int velocity_seed = (unsigned int) retrieve_item(root, default_root, "velocity_seed")->valueint;
    params->initial_velocities = 
      generate_velocities(params->temperature, velocity_seed, 
			  params->masses, params->n_dims,
			  params->n_particles);
  } else {
        params->initial_velocities = 
	  load_matrix(item->valuestring, params->n_particles, params->n_dims, 0);
  }

  //prepare output files
  item = cJSON_GetObjectItem(root, "positions_log_file");
  if(item)
    params->positions_file = fopen(item->valuestring, "w");


  item = cJSON_GetObjectItem(root, "velocities_log_file");
  if(item)
    params->positions_file = fopen(item->valuestring, "w");

  item = cJSON_GetObjectItem(root, "forces_log_file");
  if(item)
    params->positions_file = fopen(item->valuestring, "w");
    

  free(data);
  cJSON_Delete(root);
  cJSON_Delete(default_root);

  return params;
  
}

Run_Params* read_parameters_old(FILE* params_file, const Run_Params* default_params) {


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
    if(params_file == NULL)
      n = fscanf(stdin,"%s %s\n", pstrings[i][KEY], pstrings[i][VALUE]);
    else
      n = fscanf(params_file,"%s %s\n", pstrings[i][KEY], pstrings[i][VALUE]);
#ifdef DEBUG
    printf("Read [%s] and [%s]\n", pstrings[i][KEY], pstrings[i][VALUE]);
#endif
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

  
  temp_uint = process_uint(pstrings, "com_remove_period", &success_cur);
  if(success_cur)
    params->com_remove_period = temp_uint;
  else if(default_params == NULL) {
    fprintf(stderr, "Warning: assuming default com_remove_period = 1000\n");
    params->com_remove_period = 1000;
  }

  temp_uint = process_uint(pstrings, "n_dims", &success_cur);
  if(success_cur)
    params->n_dims = temp_uint;
  else if(default_params == NULL) {
    fprintf(stderr, "Could not read n_dims\n");
    exit(1);
  }

  //Now we know dimensions, get box size
  char box_size_string[12];
  params->box_size = (double*) malloc(sizeof(double) * params->n_dims);
  for(i = 0; i < params->n_dims; i++) {
    sprintf(box_size_string, "box_%d_size", i+1);
    temp_double = process_double(pstrings, box_size_string, &success_cur);
    if(success_cur)
      params->box_size[i] = temp_double;
    else if(default_params == NULL) {

#ifndef NO_PBC
      fprintf(stderr, "Could not read box_%d_size\n", i);
      exit(1);
#endif

    }
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

  //Neighbor list parameters
#ifdef NLIST
  double rcut = process_double(pstrings, "rcut", &success_cur);
  if(!success_cur) {
    fprintf(stderr, "Error: Must list rcut (cut-off)\n");
    exit(1);
  }
  double skin = process_double(pstrings, "skin", &success_cur);
  if(!success_cur) {
    skin = 0.2 * rcut;
    fprintf(stderr, "Warning: Assuming skin %g\n", skin);
  }
  Nlist_Parameters* nlist = build_nlist_params(params->n_dims, params->n_particles, params->box_size, skin, rcut);
#endif //NLIST


  //thermostat parameters
  #ifdef ANDERSON
  unsigned int seed = process_uint(pstrings, "thermostat_seed", &success_cur);
  if(!success_cur) {
    seed = 1523;
    fprintf(stderr, "Warning: Using default seed for thermostat\n");
  }
  double collision_freq = process_double(pstrings, "anderson_nu", &success_cur);
  if(!success_cur) {
    collision_freq = 10;
  }
  
  params->thermostat_parameters = build_anderson(seed, collision_freq);
  #endif

  //Bussi thermostat parameters 
  #ifdef BUSSI
  unsigned int seed = process_uint(pstrings, "thermostat_seed", &success_cur);
  if(!success_cur) {
    seed = 1523;
    fprintf(stderr, "Warning: Using default seed for thermostat\n");
  }
  double taut = process_double(pstrings, "bussi_taut", &success_cur);
  if(!success_cur) {
    taut = 0;
  }
  
  params->thermostat_parameters = build_bussi(seed, taut);
  #endif  
  

  //FF parameters
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
  params->force_parameters = build_lj(epsilon, sigma, nlist);
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
#ifdef DEBUG
      printf("Read in [%s] for [%s]\n", result, key);
#endif
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
#ifdef DEBUG
      printf("Read in [%d] for [%s]\n", result, key);
#endif
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
#ifdef DEBUG
      printf("Read in [%6f] for [%s]\n", result, key);
#endif
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
    kenergy += 0.5 * etemp * masses[i];
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
    for(j = n_dims; j < 3; j++) {
      fprintf(file, "%12g ", 0.); 
    }
      fprintf(file, "\n");
  }

  fflush(file);
  
}

void log_array(FILE* file, double* array, unsigned n_cols, unsigned n_rows, bool do_sum) {
  
  if(file == NULL)
    return;

  unsigned int i,j;
  double sum = 0;

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


double remove_com(double* velocities, double* masses, unsigned int n_dims, unsigned int n_particles) {

  unsigned int i, k;
  double com[n_dims];
  double mass_sum = 0;
  double com_mag = 0;
  
  //zero 
  for(i = 0 ; i < n_dims; i++)
    com[i] = 0;
  //calculate COM in 2 parts, sum ( mass) and sum(momentum)
  for(i = 0; i < n_particles; i++) {
    mass_sum += masses[i];
    for(k = 0; k < n_dims; k++) {
      com[k] += velocities[i * n_dims + k] / masses[i];
    }
  }

  //turn into COM
  for(i = 0; i < n_dims; i++) {
    com[i] /= mass_sum;
    com_mag += com[i] * com[i];
  }

  //remove COM motion
  for(i = 0; i < n_particles; i++) {
    for(k = 0; k < n_dims; k++) {
      velocities[i * n_dims + k] -= com[k];
    }
  }
  
  return sqrt(com_mag);
}

void free_run_params(Run_Params* params) {

#ifdef THERMOSTAT
  free_thermostat(params->thermostat_parameters);
#endif
  
  free_forces(params->force_parameters);


  free(params->initial_positions);
  free(params->initial_velocities);
  free(params->masses);
  free(params->box_size);

  free(params);

}
