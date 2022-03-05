#include "util.h"
#include "force.h"
#include "group.h"

#define PARAM_FILE_BUFFER 4096

static const char *
    default_json = " { \"com_remove_period\" : 1000, \"skin\" : 0, \
    \"seed\" : 1523, \"anderson_nu\" : 10.0, \"pressure\": 0,\
    \"temperature\": 0, \"final_positions\": \"final_positions.xyz\",\
    \"harmonic_constant\" : 1.0, \"lj_epsilon\" : 1.0, \"lj_sigma\" : 1.0, \
    \"position_log_period\" : 0, \"velocity_log_period\" : 0,\
    \"box_update_period\": 0, \"force_type\": null,\
     \"force_log_period\" : 0, \"n_images\": 1, \"langevin_gamma\": 0.1} ";

#ifdef DEBUG
const char *elements[] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti",
                          "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
                          "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
                          "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
                          "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg",
                          "Bh", "Hs", "Mt", "Ds", "Rg"};
const unsigned int n_elements = 100;
#else
const char *elements[] = {"H"};
const unsigned int n_elements = 1;
#endif

void load_json_matrix(cJSON *item, double *mat, unsigned int size, const char *message);

double *generate_velocities(double temperature, gsl_rng *rng, double *masses, unsigned int n_dims, unsigned int n_particles)
{

  double *velocities = (double *)malloc(sizeof(double) * n_dims * n_particles);

  unsigned int i, j;

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

  return (velocities);
}

void load_json(char *filename, char **data)
{

  FILE *f;

  if (!filename)
  { // stdin
    fprintf(stdout, "Assuming you'll pass parameters via stdin. Waiting...\n");
    f = stdin;
    // create buffer
    int buffer_size = PARAM_FILE_BUFFER;
    char *buffer = (char *)malloc(buffer_size);
    // read in up to buffer size
    int len = fread(buffer, sizeof(char), PARAM_FILE_BUFFER, f);
    // while we fill the buffer, keep reading
    while (len % buffer_size == 0)
    {
      buffer_size += PARAM_FILE_BUFFER;
      buffer = (char *)realloc(buffer, buffer_size);
      len += fread((void *)(buffer + len), sizeof(char), PARAM_FILE_BUFFER, f);
    }

    if (ferror(f))
    {
      perror("Error reading input stream");
      exit(1);
    }

    // trim to size
    buffer = (char *)realloc(buffer, len);
    *data = buffer;
  }
  else
  { // file
    // get file length
    f = fopen(filename, "rb");
    if (!f)
    {
      fprintf(stderr, "Failed to read file: %s\n", filename);
      exit(1);
    }
    fseek(f, 0, SEEK_END);
    long len = ftell(f);
    // reset file
    fseek(f, 0, SEEK_SET);
    // load data
    *data = (char *)malloc(len + 1);
    const size_t read_size = fread(*data, 1, len, f);
    fclose(f);
  }
}

cJSON *retrieve_item(cJSON *root, cJSON *default_root, const char *item_name)
{
  cJSON *item = cJSON_GetObjectItem(root, item_name);
  if (!item)
  { // check if it's in the default parameter list
    item = cJSON_GetObjectItem(default_root, item_name);
    if (item)
    {
      printf("Warning: assuming default value for %s = ", item_name);
      switch (item->type)
      {
      case cJSON_False:
        printf("false\n");
        break;
      case cJSON_True:
        printf("true\n");
        break;
      case cJSON_Number:
        printf("%g\n", item->valuedouble);
        break;
      case cJSON_String:
        printf("%s\n", item->valuestring);
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
  SCALAR *sdata;
  group_t *group;
  load_json(file_name, &data);
  cJSON *root = cJSON_Parse(data);
  cJSON *default_root = cJSON_Parse(default_json);
  cJSON *item, *item2;

  if (!root)
  {
    fprintf(stderr, "JSON Error before: [%s]\n", cJSON_GetErrorPtr());
    exit(1);
  }

  run_params_t *params = (run_params_t *)malloc(sizeof(run_params_t));
  params->rng = NULL;

  // get single value parameters
  params->steps = (unsigned int)retrieve_item(root, default_root, "steps")->valueint;
  params->com_remove_period = (unsigned int)retrieve_item(root, default_root, "com_remove_period")->valueint;
  params->time_step = retrieve_item(root, default_root, "time_step")->valuedouble;

  params->temperature = retrieve_item(root, default_root, "temperature")->valuedouble;

  params->n_particles = (unsigned int)retrieve_item(root, default_root, "n_particles")->valueint;

  params->print_period = (unsigned int)retrieve_item(root, default_root, "print_period")->valueint;

  params->position_log_period = (unsigned int)retrieve_item(root, default_root, "position_log_period")->valueint;
  params->position_log_period = (params->position_log_period == 0 ? params->print_period : params->position_log_period);

  params->velocity_log_period = (unsigned int)retrieve_item(root, default_root, "velocity_log_period")->valueint;
  params->velocity_log_period = (params->velocity_log_period == 0 ? params->print_period : params->velocity_log_period);

  params->force_log_period = (unsigned int)retrieve_item(root, default_root, "force_log_period")->valueint;
  params->force_log_period = (params->force_log_period == 0 ? params->print_period : params->force_log_period);

  unsigned int seed = (unsigned int)retrieve_item(root, default_root, "seed")->valueint;

  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, seed);
  params->rng = rng;

  // cell vectors
  sdata = (SCALAR *)calloc(N_DIMS * N_DIMS, sizeof(SCALAR));
  unsigned int i;
  item = cJSON_GetObjectItem(root, "cell");
  if (item)
  {
    i = 0;

    for (item = item->child; item != NULL; item = item->next)
    {
      if (i == N_DIMS)
      {
        // maybe it's default with no cell ?
        if (sdata[0] == 0)
          break;
        fprintf(stderr, "Error: Number of cell dimensions not equal to simulation dimension\n");
        exit(1);
      }
      sdata[i] = item->valuedouble;
      i++;
    }
    if (i != N_DIMS && i != N_DIMS * N_DIMS)
      fprintf(stderr, "Not enough cel dims set\n");
    // convert if was sizes, not basis vectors
    if (i == N_DIMS)
    {
      // evil wrapping
      for (i = N_DIMS - 1; i < N_DIMS; i--)
        sdata[i * N_DIMS + i] = sdata[i];
      for (i = 1; i < N_DIMS; i++)
        sdata[i] = 0;
    }
  }

  // thermostats
  params->thermostat_parameters = NULL;
  if (params->temperature)
  {
    item = cJSON_GetObjectItem(root, "thermostat");
    if (item && item->valuestring)
    {
      const char *thermostat = item->valuestring;

      if (!strcmp(thermostat, "anderson"))
      {
        double collision_freq = retrieve_item(root, default_root, "anderson_nu")->valuedouble;
        params->thermostat_parameters = build_anderson(collision_freq, params->rng);
      }
      else if (!strcmp(thermostat, "bussi"))
      {
        double taut = retrieve_item(root, default_root, "bussi_taut")->valuedouble;
        params->thermostat_parameters = build_bussi(taut, params->rng);
      }
      else if (!strcmp(thermostat, "baoab"))
      {
        double taut = retrieve_item(root, default_root, "langevin_gamma")->valuedouble;
        params->thermostat_parameters = build_baoab(taut, params->rng);
      }
      else
      {
        fprintf(stderr, "Could not understand thermostat type %s\n", thermostat);
        exit(1);
      }
    }
  }

  item = cJSON_GetObjectItem(root, "masses_file");

  // load input files
  char *positions_file = retrieve_item(root, default_root, "start_positions")->valuestring;
  params->initial_positions = load_matrix(positions_file, params->n_particles, N_DIMS, 0);
  params->scaled_positions = (double *)malloc(sizeof(double) * params->n_particles * N_DIMS);

  item = cJSON_GetObjectItem(root, "masses_file");
  if (item)
  {
    char *masses_file = item->valuestring;
    params->masses = load_matrix(masses_file, params->n_particles, 1, 0);
  }
  else
  {
    printf("Warning: Setting masses to 1.0\n");
    params->masses = (double *)malloc(sizeof(double) * params->n_particles);
    for (unsigned int i = 0; i < params->n_particles; i++)
      params->masses[i] = 1.0;
  }

  item = cJSON_GetObjectItem(root, "start_velocities");
  if (!item)
  {
    item = cJSON_GetObjectItem(root, "start_temperature");
    if (item)
      params->initial_velocities =
          generate_velocities(item->valuedouble, params->rng,
                              params->masses, N_DIMS,
                              params->n_particles);
    else
      params->initial_velocities = generate_velocities(params->temperature, params->rng,
                                                       params->masses, N_DIMS,
                                                       params->n_particles);
  }
  else
  {
    params->initial_velocities =
        load_matrix(item->valuestring, params->n_particles, N_DIMS, 0);
  }

  // read in images - first check for array
  item = cJSON_GetObjectItem(root, "images");
  unsigned int images[N_DIMS];
  if (item)
  {
    if (item->type != cJSON_Array)
    {
      fprintf(stderr, "Error: images must be an array\n");
      exit(1);
    }
    item = item->child;
    for (unsigned int i = 0; i < N_DIMS && item != NULL; i++, item = item->next)
      images[i] = item->valueint;
  }
  else
  {
    // default to n_images
    for (unsigned int i = 0; i < N_DIMS; i++)
      images[i] = retrieve_item(root, default_root, "n_images")->valueint;
  }

  // group - partition to ghost too
  // also will finish making box here, which was deferred
  group = NULL;
  item = cJSON_GetObjectItem(root, "group");
  params->dof = (params->n_particles - 1) * N_DIMS;
  if (item)
  {
    group = load_group(item->valuestring);
    group->n_gparticles = params->n_particles;
    params->n_cell_particles = group->size * params->n_particles;
    // wyckoffs are a list dictionaries - load into linked list
    item = cJSON_GetObjectItem(root, "wyckoffs");
    if (item)
    {
      group_t *next;
      unsigned int group_size = group->size;
      unsigned int w_particles = 0;
      params->n_cell_particles = 0; // reset since we are taking particles for wyckoffs
      next = group;
      for (item = item->child; item != NULL; item = item->next)
      {
        next->next = load_group(cJSON_GetObjectItem(item, "group")->valuestring);
        next = next->next;
        group_size += next->size;

        item2 = cJSON_GetObjectItem(item, "n_particles");
        if (!item2)
        {
          fprintf(stderr, "Error: n_particles must be set for each wyckoff\n");
          exit(1);
        }
        next->n_gparticles = item2->valueint;
        w_particles += next->n_gparticles;
      }
      if (w_particles > params->n_particles)
      {
        fprintf(stderr, "Error: n_particles in wyckoffs exceeds n_particles\n");
        exit(1);
      }
      group->n_gparticles = params->n_particles - w_particles;

      // for simplicity, set total size at each group
      // Also, count number of ghost particles (those not in asymmetric unit)
      for (next = group; next != NULL; next = next->next)
      {
        next->total_size = group_size;
        params->dof += next->n_gparticles * next->dof;
        params->n_cell_particles += next->n_gparticles * next->size;
#ifdef DEBUG
        printf("Loaded group %s with %d particles and %d members\n", next->name, next->n_gparticles, next->size);
#endif
      }
      // remove translation dof
      params->dof -= N_DIMS;
    }
    params->box = make_box(sdata, group, images);
    sdata = NULL;
    // expand ghost particles to include tilings
    // we subtract at the end to avoid counting the asymmetric unit particles
    params->n_ghost_particles = params->n_cell_particles * (params->box->n_tilings + 1) - params->n_particles;
#ifdef DEBUG
    printf("Duplicating %d particles into %d real particles and %d ghost for group with %d elements and %d tilings. Each cell has %d particles.\n",
           params->n_particles, params->n_particles, params->n_ghost_particles, params->box->group->size, params->box->n_tilings, params->n_cell_particles);
    printf("Computed %d degrees of freedom\n", params->dof);
#endif

    // make new longer list
    sdata = (SCALAR *)malloc(sizeof(SCALAR) * N_DIMS * (params->n_ghost_particles + params->n_particles));

    // before we're done, use it to unscale
    for (i = 0; i < params->n_particles; i++)
      unscale_coords(&sdata[i * N_DIMS],
                     &params->initial_positions[i * N_DIMS], params->box);
    // ok now done
    free(params->initial_positions);
    params->initial_positions = sdata;
    sdata = NULL;

    // now need velocities and msses
    // make new longer list
    sdata = (SCALAR *)malloc(sizeof(SCALAR) * N_DIMS * (params->n_ghost_particles + params->n_particles));
    memcpy(sdata, params->initial_velocities, params->n_particles * sizeof(SCALAR) * N_DIMS);
    // before we're done, use it to unscale
    // ok now done
    free(params->initial_velocities);
    params->initial_velocities = sdata;
    sdata = NULL;

    // TODO: This code is not correct.  It consider Wyckoff groups that are different size
    sdata = (SCALAR *)calloc(params->n_cell_particles, sizeof(SCALAR));
    for (i = 0; i + params->n_particles <= params->n_cell_particles; i += params->n_particles)
      memcpy(&sdata[i], params->masses, params->n_particles * sizeof(SCALAR));
    // before we're done, use it to unscale
    // ok now done
    free(params->masses);
    params->masses = sdata;
    sdata = NULL;
  }
  else
  {
    params->box = make_box(sdata, NULL, images);
    sdata = NULL;
  }

  params->box_update_period = (unsigned int)retrieve_item(root, default_root, "box_update_period")->valueint;
  if (params->box_update_period)
  {
    params->pressure = (double)retrieve_item(root, default_root, "pressure")->valuedouble;
  }

  // forces
  const char *force_type = retrieve_item(root, default_root, "force_type")->valuestring;
  if (!force_type)
  {
    params->force_parameters = NULL;
  }
  else if (!strcmp(force_type, "harmonic"))
  {
    double k = retrieve_item(root, default_root, "harmonic_constant")->valuedouble;
    params->force_parameters = build_harmonic(k);
  }
  else if (!strcmp(force_type, "nlj"))
  {
    double epsilon = retrieve_item(root, default_root, "lj_epsilon")->valuedouble;
    double sigma = retrieve_item(root, default_root, "lj_sigma")->valuedouble;
    // make nlist
    nlist_parameters_t *nlist = NULL;
    double rcut = sigma * 3;
    double skin = 0.2 * rcut;
    item = cJSON_GetObjectItem(root, "rcut");
    if (item)
    {
      rcut = item->valuedouble;
      skin = retrieve_item(root, default_root, "skin")->valuedouble;
      if (skin == 0)
      {
        skin = 0.2 * rcut;
        printf("Warning: Assuming skin = %g\n", skin);
      }
    }

    nlist = build_nlist_params(N_DIMS, params->n_particles, params->n_ghost_particles,
                               params->box->box_size, skin, rcut);
    params->force_parameters = build_nlj(epsilon, sigma, nlist);
  }
  else if (!strcmp(force_type, "lj"))
  {
    double epsilon = retrieve_item(root, default_root, "lj_epsilon")->valuedouble;
    double sigma = retrieve_item(root, default_root, "lj_sigma")->valuedouble;
    params->force_parameters = build_lj(epsilon, sigma);
  }
  else if (!strcmp(force_type, "soft"))
  {
    params->force_parameters = build_soft();
  }
  else if (!strcmp(force_type, "gravity"))
  {
    params->force_parameters = build_gravity();
  }
  else
  {
    fprintf(stderr, "Could not understand force type %s\n", force_type);
    exit(1);
  }

  // final frame
  item = retrieve_item(root, default_root, "final_positions");
  params->final_positions_file = fopen(item->valuestring, "w");

  // prepare output files
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

group_t *load_group(char *filename)
{
  char *data;
  unsigned g_dims = N_DIMS + 1;
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

  item = cJSON_GetObjectItem(root, "dof");
  if (!item)
  {
    fprintf(stderr, "Malformed group JSON - must have dof\n");
    exit(1);
  }
  group->dof = item->valueint;

  unsigned int size = 0;
  g_t *tmp, *members = NULL;
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
    // load member
    load_json_matrix(json_members, members[size - 1].g, g_dims * g_dims, "g matrix");
  }
  group->members = members;
  group->size = size;
  group->total_size = size;

  load_json_matrix(cJSON_GetObjectItem(root, "projector"), group->projector, N_DIMS * N_DIMS * N_DIMS * N_DIMS, "projector");
  group->next = NULL;

  free(data);

  return group;
}

void load_json_matrix(cJSON *item, double *mat, unsigned int size, const char *message)
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
      fprintf(stderr, "Too many values in %s. Only found %d out of %d\n", message, i, size);
      exit(1);
    }
    mat[i++] = item->valuedouble;
  }
  if (i != size)
  {
    fprintf(stderr, "Too few values in %s. Only found %d out of %d\n", message, i, size);
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
      etemp += velocities[i * n_dims + j] * velocities[i * n_dims + j];
    kenergy += 0.5 * etemp * masses[i];
  }

  return (kenergy);
}

void log_xyz(FILE *file, double *array, char *frame_string,
             unsigned n_dims, unsigned n_particles,
             unsigned int total, int location)
{

  if (file == NULL)
    return;

  unsigned int i, j;
  if (location == 0)
    fprintf(file, "%d\n%s\n", total, frame_string);

  for (i = 0; i < n_particles; i++)
  {
    fprintf(file, "%s ", elements[i % n_elements]);
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
  if (location == 2)
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
    fprintf(stderr, "file: %s", filename);
    exit(1);
  }

  return NULL;
}

double remove_com(double *data, double *masses, unsigned int n_dims, unsigned int n_particles)
{

  unsigned int i, k;
  double com[n_dims];
  double mass_sum = 0;
  double com_mag = 0;

  // zero
  for (i = 0; i < n_dims; i++)
    com[i] = 0;
  // calculate COM in 2 parts, sum ( mass) and sum(momentum)
  for (i = 0; i < n_particles; i++)
  {
    mass_sum += masses[i];
    for (k = 0; k < n_dims; k++)
    {
      com[k] += data[i * n_dims + k] / masses[i];
    }
  }

  // turn into COM
  for (i = 0; i < n_dims; i++)
  {
    com[i] /= mass_sum;
    com_mag += com[i] * com[i];
  }

  // remove COM motion
  for (i = 0; i < n_particles; i++)
  {
    for (k = 0; k < n_dims; k++)
    {
      data[i * n_dims + k] -= com[k];
    }
  }

  return sqrt(com_mag);
}

void free_run_params(run_params_t *params)
{

  if (params->thermostat_parameters)
    params->thermostat_parameters->free(params->thermostat_parameters);
  if (params->force_parameters)
    params->force_parameters->free(params->force_parameters);

  free(params->initial_positions);
  free(params->initial_velocities);
  free(params->masses);
  free_box(params->box);
  if (params->positions_file)
    fclose(params->positions_file);
  if (params->velocities_file)
    fclose(params->velocities_file);
  if (params->forces_file)
    fclose(params->forces_file);

  gsl_rng_free(params->rng);

  free(params);
}
