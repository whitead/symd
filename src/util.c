#include "util.h"
#include "force.h"
#include "group.h"

#define PARAM_FILE_BUFFER 1024

static const char *
    default_json = " { \"com_remove_period\" : 1000, \"skin\" : 0, \
    \"thermostat_seed\" : 1523, \"anderson_nu\" : 10.0,\
    \"temperature\": 0, \"rcut\": 0, \"final_positions\": \"final_positions.xyz\",\
    \"harmonic_constant\" : 1.0, \"lj_epsilon\" : 1.0, \"lj_sigma\" : 1.0, \
    \"velocity_seed\" : 543214, \"position_log_period\" : 0, \"velocity_log_period\" : 0,\
     \"force_log_period\" : 0, \"box_size\": [0, 0, 0]} ";

void *load_json_matrix(cJSON *item, double *mat, unsigned int size, const char *message);

double *generate_velocities(double temperature, unsigned int seed, double *masses, unsigned int n_dims, unsigned int n_particles)
{

  double *velocities = (double *)malloc(sizeof(double) * n_dims * n_particles);

  unsigned int i, j;
  gsl_rng *rng;

  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, seed);

  for (i = 0; i < n_particles; i++)
  {
    for (j = 0; j < n_dims; j++)
    {
      if (temperature != 0)
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

  return (velocities);
}

void load_json(char *filename, char **data)
{

  FILE *f;

  if (!filename)
  { //stdin
    fprintf(stdout, "Assuming you'll pass parameters via stdin. Waiting...\n");
    f = stdin;
    //create buffer
    int buffer_size = PARAM_FILE_BUFFER;
    char *buffer = (char *)malloc(buffer_size);
    //read in up to buffer size
    int len = fread(buffer, sizeof(char), PARAM_FILE_BUFFER, f);
    //while we fill the buffer, keep reading
    while (len % buffer_size == 0)
    {
      buffer_size += PARAM_FILE_BUFFER;
      buffer = (char *)realloc(buffer, buffer_size);
      len += fread((void *)(buffer + len), sizeof(char), PARAM_FILE_BUFFER, f);
    }

    if (ferror(f))
    {
      perror("Error reading input stream");
    }

    //trim to size
    buffer = (char *)realloc(buffer, len);
    *data = buffer;
  }
  else
  { //file
    //get file length
    f = fopen(filename, "rb");
    fseek(f, 0, SEEK_END);
    long len = ftell(f);
    //reset file
    fseek(f, 0, SEEK_SET);
    //load data
    *data = (char *)malloc(len + 1);
    const size_t read_size = fread(*data, 1, len, f);
    fclose(f);
  }
}

cJSON *retrieve_item(cJSON *root, cJSON *default_root, const char *item_name)
{
  cJSON *item = cJSON_GetObjectItem(root, item_name);
  if (!item)
  { //check if it's in the default parameter list
    item = cJSON_GetObjectItem(default_root, item_name);
    if (item)
    {
      fprintf(stderr, "Warning: assuming default value for %s = ", item_name);
      switch (item->type)
      {
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
    else
    {
      fprintf(stderr, "Error: could not read %s and no default value\n", item_name);
      exit(1);
    }
  }

  return item;
}

run_params_t *read_parameters(char *file_name)
{

  char *data;
  load_json(file_name, &data);
  cJSON *root = cJSON_Parse(data);
  cJSON *default_root = cJSON_Parse(default_json);
  cJSON *item;

  if (!root)
  {
    fprintf(stderr, "Error before: [%s]\n", cJSON_GetErrorPtr());
    exit(1);
  }

  run_params_t *params = (run_params_t *)malloc(sizeof(run_params_t));

  //get single value parameters
  params->steps = (unsigned int)retrieve_item(root, default_root, "steps")->valueint;
  params->com_remove_period = (unsigned int)retrieve_item(root, default_root, "com_remove_period")->valueint;
  params->time_step = retrieve_item(root, default_root, "time_step")->valuedouble;

  params->temperature = retrieve_item(root, default_root, "temperature")->valuedouble;

  params->n_dims = (unsigned int)retrieve_item(root, default_root, "n_dims")->valueint;
  params->n_particles = (unsigned int)retrieve_item(root, default_root, "n_particles")->valueint;

  params->print_period = (unsigned int)retrieve_item(root, default_root, "print_period")->valueint;

  params->position_log_period = (unsigned int)retrieve_item(root, default_root, "position_log_period")->valueint;
  params->position_log_period = (params->position_log_period == 0 ? params->print_period : params->position_log_period);

  params->velocity_log_period = (unsigned int)retrieve_item(root, default_root, "velocity_log_period")->valueint;
  params->velocity_log_period = (params->velocity_log_period == 0 ? params->print_period : params->velocity_log_period);

  params->force_log_period = (unsigned int)retrieve_item(root, default_root, "force_log_period")->valueint;
  params->force_log_period = (params->force_log_period == 0 ? params->print_period : params->force_log_period);

  params->box_size = (double *)calloc(params->n_dims, sizeof(double));
  //box size
  unsigned int i;
  item = retrieve_item(root, default_root, "box_size");
  i = 0;
  for (item = item->child; item != NULL; item = item->next)
  {
    if (i == params->n_dims)
    {
      // maybe it's default with no box?
      if (params->box_size[0] == 0)
        break;
      fprintf(stderr, "Error: Number of box dimensions not equal to simulation dimension\n");
      exit(1);
    }
    params->box_size[i] = item->valuedouble;
    i++;
  }
  //thermostats
  params->thermostat_parameters = NULL;
  if (params->temperature)
  {
    item = cJSON_GetObjectItem(root, "thermostat");
    if (item)
    {
      const char *thermostat = item->valuestring;

      if (!strcmp(thermostat, "anderson"))
      {
        unsigned int seed = (unsigned int)retrieve_item(root, default_root, "thermostat_seed")->valueint;
        double collision_freq = retrieve_item(root, default_root, "anderson_nu")->valuedouble;
        params->thermostat_parameters = build_anderson(seed, collision_freq);
      }
      else if (!strcmp(thermostat, "bussi"))
      {
        unsigned int seed = (unsigned int)retrieve_item(root, default_root, "thermostat_seed")->valueint;
        double taut = retrieve_item(root, default_root, "bussi_taut")->valuedouble;
        params->thermostat_parameters = build_bussi(seed, taut);
      }
      else
      {
        fprintf(stderr, "Could not understand thermostat type %s\n", thermostat);
        exit(1);
      }
    }
  }

  //load input files
  char *positions_file = retrieve_item(root, default_root, "start_positions")->valuestring;
  params->initial_positions = load_matrix(positions_file, params->n_particles, params->n_dims, 0);

  item = cJSON_GetObjectItem(root, "masses_file");
  if (item)
  {
    char *masses_file = item->valuestring;
    params->masses = load_matrix(masses_file, params->n_particles, 1, 0);
  }
  else
  {
    fprintf(stderr, "Warning: Setting masses to 1.0\n");
    params->masses = (double *)malloc(sizeof(double) * params->n_particles);
    for (unsigned int i = 0; i < params->n_particles; i++)
      params->masses[i] = 1.0;
  }

  item = cJSON_GetObjectItem(root, "start_velocities");
  if (!item)
  {
    unsigned int velocity_seed = (unsigned int)retrieve_item(root, default_root, "velocity_seed")->valueint;
    params->initial_velocities =
        generate_velocities(params->temperature, velocity_seed,
                            params->masses, params->n_dims,
                            params->n_particles);
  }
  else
  {
    params->initial_velocities =
        load_matrix(item->valuestring, params->n_particles, params->n_dims, 0);
  }

  // group - partition to ghost too
  item = cJSON_GetObjectItem(root, "group");
  if (item)
  {
    params->group = load_group(item->valuestring, params->n_dims);
    params->n_particles = params->n_particles / params->group->size;
    params->n_ghost_particles = params->n_particles * (params->group->size - 1);
    // remove velocties on ghost particles
    for (unsigned int i = params->n_particles; i < params->n_particles + params->n_ghost_particles; i++)
      for (unsigned int j = 0; j < params->n_dims; j++)
        params->initial_velocities[i * params->n_dims + j] = 0;
  }
  else
  {
    params->group = NULL;
    params->n_ghost_particles = 0;
  }

  // make nlist
  nlist_parameters_t *nlist = NULL;
  item = retrieve_item(root, default_root, "rcut");
  if (item)
  {
    double rcut = item->valuedouble;
    double skin = retrieve_item(root, default_root, "skin")->valuedouble;
    if (skin == 0)
    {
      skin = 0.2 * rcut;
      fprintf(stderr, "Warning: Assuming skin = %g\n", skin);
    }
    nlist = build_nlist_params(params->n_dims, params->n_particles, params->n_ghost_particles,
                               params->box_size, skin, rcut);
  }

  // forces
  const char *force_type = retrieve_item(root, default_root, "force_type")->valuestring;
  if (!strcmp(force_type, "harmonic"))
  {
    double k = retrieve_item(root, default_root, "harmonic_constant")->valuedouble;
    params->force_parameters = build_harmonic(k);
  }
  else if (!strcmp(force_type, "lj"))
  {
    double epsilon = retrieve_item(root, default_root, "lj_epsilon")->valuedouble;
    double sigma = retrieve_item(root, default_root, "lj_sigma")->valuedouble;
    params->force_parameters = build_lj(epsilon, sigma, nlist);
  }
  else if (!strcmp(force_type, "soft"))
  {
    params->force_parameters = build_soft();
  }
  else
  {
    fprintf(stderr, "Could not understand force type %s\n", force_type);
    exit(1);
  }

  // final frame
  item = retrieve_item(root, default_root, "final_positions");
  params->final_positions_file = fopen(item->valuestring, "w");

  //prepare output files
  item = cJSON_GetObjectItem(root, "positions_log_file");
  if (item)
    params->positions_file = fopen(item->valuestring, "w");
  else
    params->positions_file = NULL;

  item = cJSON_GetObjectItem(root, "velocities_log_file");
  if (item)
    params->velocities_file = fopen(item->valuestring, "w");
  else
    params->velocities_file = NULL;

  item = cJSON_GetObjectItem(root, "forces_log_file");
  if (item)
    params->forces_file = fopen(item->valuestring, "w");
  else
    params->forces_file = NULL;

  free(data);
  cJSON_Delete(root);
  cJSON_Delete(default_root);

  return params;
}

group_t *load_group(char *filename, unsigned int n_dims)
{
  char *data;
  load_json(filename, &data);
  cJSON *root = cJSON_Parse(data);
  cJSON *item;

  if (!root)
  {
    fprintf(stderr, "Error in group JSON before: [%s]\n", cJSON_GetErrorPtr());
    exit(1);
  }

  item = cJSON_GetObjectItem(root, "name");
  cJSON *json_members = cJSON_GetObjectItem(root, "members");
  if (!json_members || !item)
  {
    fprintf(stderr, "Malformed group JSON - must have members and name\n");
    exit(1);
  }

  group_t *group = (group_t *)malloc(sizeof(group_t));
  group->name = item->valuestring;
  unsigned int i, j, size = 0;
  g_t *tmp, *members = NULL;
  double *g_mat, *i_mat;
  // iterate over group members

  for (json_members = json_members->child; json_members != NULL; json_members = json_members->next)
  {

    // grow
    tmp = (g_t *)malloc(sizeof(g_t) * ++size);
    if (members)
    {
      memcpy(tmp, members, sizeof(g_t) * (size - 1));
      free(members);
    }
    members = tmp;

    g_mat = (double *)malloc(sizeof(double) * n_dims * n_dims);
    i_mat = (double *)malloc(sizeof(double) * n_dims * n_dims);
    load_json_matrix(cJSON_GetObjectItem(json_members, "g"), g_mat, n_dims * n_dims, "g matrix");
    load_json_matrix(cJSON_GetObjectItem(json_members, "i"), i_mat, n_dims * n_dims, "i matrix");
    g_t g = {.g = g_mat, .i = i_mat};

    // add new member
    memcpy(&members[size - 1], &g, sizeof(g_t));
  }
  group->members = members;
  group->size = size;

  free(data);

  return group;
}

void *load_json_matrix(cJSON *item, double *mat, unsigned int size, const char *message)
{
  if (!item)
  {
    fprintf(stderr, "Malformed group JSON - could not parse %s\n", message);
    exit(1);
  }
  unsigned int i = 0;
  for (item = item->child; item != NULL; item = item->next)
  {
    if (i == size)
    {
      fprintf(stderr, "Too many values in %s\n", message);
      exit(1);
    }
    mat[i++] = item->valuedouble;
  }
  if (i != size)
  {
    fprintf(stderr, "Too few values in %s\n", message);
    exit(1);
  }
}

double calculate_kenergy(double *velocities, double *masses, unsigned int n_dims, unsigned int n_particles)
{

  unsigned int i, j;
  double kenergy = 0;
  double etemp;
  for (i = 0; i < n_particles; i++)
  {
    etemp = 0;
    for (j = 0; j < n_dims; j++)
    {
      etemp += velocities[i * n_dims + j] * velocities[i * n_dims + j];
    }
    kenergy += 0.5 * etemp * masses[i];
  }

  return (kenergy);
}

void log_xyz(FILE *file, double *array, char *frame_string, unsigned n_dims, unsigned n_particles)
{

  if (file == NULL)
    return;

  unsigned int i, j;

  fprintf(file, "%d\n%s\n", n_particles, frame_string);

  for (i = 0; i < n_particles; i++)
  {
    fprintf(file, "Ar ");
    for (j = 0; j < n_dims; j++)
    {
      fprintf(file, "%12g ", array[i * n_dims + j]);
    }
    for (j = n_dims; j < 3; j++)
    {
      fprintf(file, "%12g ", 0.);
    }
    fprintf(file, "\n");
  }

  fflush(file);
}

void log_array(FILE *file, double *array, unsigned n_cols, unsigned n_rows, bool do_sum)
{

  if (file == NULL)
    return;

  unsigned int i, j;
  double sum = 0;

  for (i = 0; i < n_rows; i++)
  {
    if (do_sum)
      sum = 0;
    for (j = 0; j < n_cols; j++)
    {
      if (do_sum)
        sum += array[i * n_cols + j] * array[i * n_cols + j];
      fprintf(file, "%12g ", array[i * n_cols + j]);
    }
    if (do_sum)
      fprintf(file, "%12g\n", sqrt(sum));
    else
      fprintf(file, "\n");
  }

  fflush(file);
}

double *load_matrix(char *filename, unsigned int nrow, unsigned int ncol, unsigned int skip)
{

  FILE *mfile;
  mfile = fopen(filename, "r");
  if (mfile != NULL)
  {

    int c;
    unsigned int i, j;

    if (skip > 0)
    {
      for (i = 0; i < skip; i++)
      {
        while ((c = getc(mfile)) != '\n' && c != EOF)
          ;
      }
    }

    double *matrix = (double *)malloc(sizeof(double) * nrow * ncol);
    if (matrix == NULL)
    {
      perror("Out of memory\n");
      exit(1);
    }

    for (i = 0; i < nrow; i++)
    {
      for (j = 0; j < ncol; j++)
      {
        if (fscanf(mfile, "%lg", &(matrix[i * ncol + j])) == 0)
        {
          fprintf(stderr, "Incorrect number of rows or columns"
                          "at i = %d, and j=%d, nrow=%d, ncol=%d\n",
                  i, j, nrow, ncol);
          exit(1);
        }
      }
    }

    fclose(mfile);
    return matrix;
  }
  else
  {
    perror("Could not open file\n");
  }

  return NULL;
}

double remove_com(double *velocities, double *masses, unsigned int n_dims, unsigned int n_particles)
{

  unsigned int i, k;
  double com[n_dims];
  double mass_sum = 0;
  double com_mag = 0;

  //zero
  for (i = 0; i < n_dims; i++)
    com[i] = 0;
  //calculate COM in 2 parts, sum ( mass) and sum(momentum)
  for (i = 0; i < n_particles; i++)
  {
    mass_sum += masses[i];
    for (k = 0; k < n_dims; k++)
    {
      com[k] += velocities[i * n_dims + k] / masses[i];
    }
  }

  //turn into COM
  for (i = 0; i < n_dims; i++)
  {
    com[i] /= mass_sum;
    com_mag += com[i] * com[i];
  }

  //remove COM motion
  for (i = 0; i < n_particles; i++)
  {
    for (k = 0; k < n_dims; k++)
    {
      velocities[i * n_dims + k] -= com[k];
    }
  }

  return sqrt(com_mag);
}

void free_run_params(run_params_t *params)
{

  if (params->thermostat_parameters)
    params->thermostat_parameters->free(params->thermostat_parameters);
  params->force_parameters->free(params->force_parameters);

  free(params->initial_positions);
  free(params->initial_velocities);
  free(params->masses);
  free(params->box_size);
  if (params->group)
    free_group(params->group);

  if (params->positions_file)
    fclose(params->positions_file);
  if (params->velocities_file)
    fclose(params->velocities_file);
  if (params->forces_file)
    fclose(params->forces_file);

  free(params);
}
