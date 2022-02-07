#include "nlist.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//TODO: Should be 2.0 in principle? - Check this
#define BOX_M 2.0

/*
 * Grow and insert a value as needed
 */
unsigned int insert_grow(unsigned int index, unsigned int value, unsigned int **array, unsigned int length);

/*
 * build the neighbor list
 */
void build_list(double *positions, double *box_size, unsigned int n_dims, unsigned int n_particles, unsigned int n_ghost_particles, nlist_parameters_t *nlist);

/*
 * Precompute cells used for constructing neighbor lists.
 * Note: this is not ready to be called multiple times. If there
 * ever is barostatting, this needs work.
 */
void build_cells(double *box_size, unsigned int n_dims,
                 unsigned int n_particles, unsigned int n_ghost_particles,
                 nlist_parameters_t *nlist);

nlist_parameters_t *build_nlist_params(unsigned int n_dims, unsigned int n_particles, unsigned int n_ghost_particles, double *box_size, double skin, double rcut)
{

  nlist_parameters_t init = {.rcut = rcut * rcut, .skin = skin * skin, .skin_rcut = (rcut + skin) * (rcut + skin)};
  nlist_parameters_t *nlist = (nlist_parameters_t *)malloc(sizeof(nlist_parameters_t));
  memcpy(nlist, &init, sizeof(nlist_parameters_t));
  nlist->nlist = NULL;
  nlist->last_positions = (double *)malloc(sizeof(double) * n_particles * n_dims);
  build_cells(box_size, n_dims, n_particles, n_ghost_particles, nlist);
  return nlist;
}

void free_nlist(nlist_parameters_t *nlist)
{

  free(nlist->last_positions);
  free(nlist->nlist);
  free(nlist->nlist_count);
  free(nlist->cell_number);
  free(nlist->adjacent_cells);
  free(nlist->mapping);
  free(nlist->head);
  free(nlist->cell_list);

  free(nlist);
}

void update_nlist(double *positions,
                  double *box_size,
                  unsigned int n_dims,
                  unsigned int n_particles,
                  unsigned int n_ghost_particles,
                  nlist_parameters_t *nlist)
{

  if (nlist->nlist == NULL)
  {
    build_list(positions, box_size, n_dims, n_particles, n_ghost_particles, nlist);
  }

  double dist, temp;
  double max1, max2;
  max1 = max2 = 0;
  unsigned int i, k;

  for (i = 0; i < n_particles; i++)
  {
    dist = 0;
    for (k = 0; k < n_dims; k++)
    {
      //temp = (positions[i * n_dims + k] - nlist->last_positions[i * n_dims + k]);
      temp = min_image_dist(temp, box_size[k]);
      dist += temp * temp;
    }
    if (dist > max1)
    {
      max1 = dist;
    }
    else if (dist > max2)
    {
      max2 = dist;
    }
  }
  if (max1 + max2 > nlist->skin)
  {
#ifdef DEBUG
    printf("updating nlist due to %f + %f > %f\n", max1, max2, nlist->skin);
#endif
    //build_list(positions, box_size, n_dims, n_particles, n_ghost_particles, nlist);
  }
  //**********************
  //TODO: This line below should be removed
  //**********************
  build_list(positions, box_size, n_dims, n_particles, n_ghost_particles, nlist);
  return;
}

unsigned int gen_index_r(int *array, int *dims, unsigned int index, unsigned int cur_dim, unsigned int n_dims, int value, unsigned int output_length)
{
  int diff[3] = {0, -1, 1};
  unsigned int i;
  for (i = 0; index < output_length && i < 3; i++)
  {
    if (cur_dim < n_dims - 1)
      index = gen_index_r(array, dims, index, cur_dim + 1, n_dims, value * dims[cur_dim] + diff[i], output_length);
    else
    {
      array[index] = value * dims[cur_dim] + diff[i];
      index++;
    }
  }

  return index;
}

void build_cells(double *box_size,
                 unsigned int n_dims,
                 unsigned int n_particles,
                 unsigned int n_ghost_particles,
                 nlist_parameters_t *nlist)
{

  //rebuild the list
  unsigned int i;
  nlist->cell_number = (int *)malloc(sizeof(int) * n_dims);
  double m_box_size[n_dims];

  // make our box bigger, since our tiled ghost particles
  // will not wrap nicely.
  memcpy(m_box_size, box_size, n_dims * sizeof(double));
  for (i = 0; i < n_dims; i++)
    m_box_size[i] *= BOX_M;

  //build cell list
  nlist->cell_number_total = 1;
  double width = sqrt(nlist->skin_rcut);
  for (i = 0; i < n_dims; i++)
  {
    nlist->cell_number[i] =
        (int)(m_box_size[i] + width - 1) / width; // ceiling quotient
    //minimum
    if (nlist->cell_number[i] == 0)
      nlist->cell_number[i] = 1;
    nlist->cell_number[i] += 2; //for ghost cells
    nlist->cell_number_total *=
        nlist->cell_number[i] - 2; //don't include ghost cells in calculation
  }

  nlist->ncell_number = (pow(3, n_dims) - 1) / 2 + 1; //number of neighboring cells
  if (nlist->ncell_number > nlist->cell_number_total) //for small systems
    nlist->ncell_number = nlist->cell_number_total;

#ifdef DEBUG
  printf("cell number total: %d (", nlist->cell_number_total);
  for (i = 0; i < n_dims; i++)
    printf("%d ", nlist->cell_number[i]);
  printf(")\ncell number total: %d\n", nlist->cell_number_total);
  printf("neighbor of neighbor cells: %d\n", nlist->ncell_number);
#endif //DEBUG

  nlist->adjacent_cells = (int *)malloc(sizeof(int) * nlist->ncell_number);

  //create offset matrix
  gen_index_r(nlist->adjacent_cells, nlist->cell_number, 0, 0, n_dims, 0, nlist->ncell_number);

#ifdef DEBUG
  printf("The neighbors offset:\n");
  for (i = 0; i < nlist->ncell_number; i++)
    printf("%d\n", nlist->adjacent_cells[i]);
#endif //DEBUG

  //create mapping for PBCs
  //find how far the cells may extend past real cells with the offset matrix
  nlist->map_offset = 0;
  for (i = 0; i < n_dims; i++)
    nlist->map_offset = nlist->map_offset * nlist->cell_number[i] + 1;

  nlist->mapping =
      (unsigned int *)malloc(sizeof(unsigned int) *
                             (nlist->map_offset * 2 + nlist->cell_number_total));

  //wrapped mapping
  for (i = 0; i < nlist->map_offset; i++)
  {
    nlist->mapping[i] =
        nlist->cell_number_total - i % nlist->cell_number_total - 1;
    nlist->mapping[i + nlist->map_offset + nlist->cell_number_total] =
        i % nlist->cell_number_total;
  }

  //trivial mapping
  for (i = 0; i < nlist->cell_number_total; i++)
    nlist->mapping[i + nlist->map_offset] = i;

#ifdef DEBUG
  printf("mapping:\n");
  for (i = 0; i < nlist->map_offset * 2 + nlist->cell_number_total; i++)
    printf("%d ", nlist->mapping[i]);
  printf("\n");
#endif //DEBUG

  //remove ghost cells from cell count
  for (i = 0; i < n_dims; i++)
    nlist->cell_number[i] -= 2;

  //prepare cell list and head
  nlist->head = (int *)malloc(sizeof(int) * nlist->cell_number_total);
  nlist->cell_list = (int *)malloc(sizeof(int) * (n_particles + n_ghost_particles));
}

/*
 * Build neighbor list using cells in N-dimensions
 */
void build_list(double *positions, double *box_size, unsigned int n_dims, unsigned int n_particles, unsigned int n_ghost_particles, nlist_parameters_t *nlist)
{

  unsigned int total_n = 0;

  if (nlist->nlist == NULL)
  {
    //TODO: Maybe only need n_particles?
    unsigned int ps = n_particles + n_ghost_particles;
    //first call, need to set-up some things
    nlist->nlist_count = (unsigned int *)malloc(sizeof(unsigned int) * ps);
    nlist->nlist = (unsigned int *)malloc(sizeof(unsigned int) * ps * (ps / 2));
    nlist->nlist_length = ps * (ps / 2);
  }

  //rebuild the list
  unsigned int i, k;
  int j;
  unsigned int ncell;
  double dist, dx;
  int net_dcell;
  int icell;
  double m_box_size[n_dims];

  // make our box bigger, since our tiled ghost particles
  // will not wrap nicely.
  memcpy(m_box_size, box_size, n_dims * sizeof(double));
  for (i = 0; i < n_dims; i++)
    m_box_size[i] *= BOX_M;

  //reset cell head
  for (i = 0; i < nlist->cell_number_total; i++)
    nlist->head[i] = -1;

  //fills the list in reverse order. I didn't invent this algorithm. The person who did is a genius.
  for (i = 0; i < n_particles + n_ghost_particles; i++)
  {
    // the cell for particle i
    icell = 0;
    for (k = 0; k < n_dims; k++)
    {
      // get 1D index of cell
      icell = (int)(wrap(positions[i * n_dims + k], m_box_size[k]) / m_box_size[k]) * nlist->cell_number[k] + icell * nlist->cell_number[k];
    }
    //have cell_list point to cell for particle i
    nlist->cell_list[i] = nlist->head[icell];
    // point this cell at particle i
    // so we can loop via this
    nlist->head[icell] = i;
  }

#ifdef DEBUG
  printf("cell_list per particle:\n");
  for (i = 0; i < n_particles + n_ghost_particles; i++)
    printf("%d ", nlist->cell_list[i]);

  printf("\nhead for each cell:\n");
  for (i = 0; i < nlist->cell_number_total; i++)
    printf("%d ", nlist->head[i]);
  printf("\n");
#endif //DEBUG

  // TODO: maybe only need n_particles???
  for (i = 0; i < n_particles + n_ghost_particles; i++)
  {
    //reset count
    nlist->nlist_count[i] = 0;

    //loop over neighbor cells, with no double counting
    icell = 0;
    //find index of particle using polynomial indexing
    for (k = 0; k < n_dims; k++)
      icell = (int)(wrap(positions[i * n_dims + k], m_box_size[k]) / m_box_size[k]) * nlist->cell_number[k] + icell * nlist->cell_number[k];
    //loop over neighbor cells
    for (ncell = 0; ncell < nlist->ncell_number; ncell++)
    {

      net_dcell = icell + nlist->adjacent_cells[ncell];
      //get head of list, using mapping
      j = nlist->head[nlist->mapping[net_dcell + nlist->map_offset]];

      //-1 marks end of cell list
      while (j != -1)
      {

        if (nlist->mapping[net_dcell + nlist->map_offset] == 0 && i >= (unsigned int)j)
        {
          //don't double count pairs in self-box
          break;
        }

        dist = 0;
        for (k = 0; k < n_dims; k++)
        {
          dx = positions[i * n_dims + k] - positions[j * n_dims + k];
          dx = min_image_dist(dx, m_box_size[k]);
          dist += dx * dx;
        }

        if (dist < nlist->skin_rcut)
        {
          nlist->nlist_length = insert_grow(total_n, j, &(nlist->nlist), nlist->nlist_length);
          total_n += 1;
          nlist->nlist_count[i] += 1;
        }

        //get next cell member
        j = nlist->cell_list[j];
      }
    }
  }

#ifdef DEBUG
  printf("nlist: \n");
  unsigned int count = 0;
  for (i = 0; i < n_particles; i++)
  {
    printf("nlist_count[%d] = %d\n\t", i, nlist->nlist_count[i]);
    for (j = 0; j < nlist->nlist_count[i]; j++)
    {
      printf("%d ", nlist->nlist[count]);
      count++;
    }
    printf("\n");
  }
  printf("Total particles: %d\n", count);
#endif //DEBUG

  //Cache old positions
  //TODO: Maybe only need n_particles?
  for (i = 0; i < n_particles; i++)
    for (k = 0; k < n_dims; k++)
      nlist->last_positions[i * n_dims + k] = positions[i * n_dims + k];
}

/* An insert that grows the array if necessary
 *
 */
unsigned int insert_grow(unsigned int index, unsigned int value, unsigned int **array, unsigned int length)
{

  if (index >= length)
  {
    unsigned int new_length = (unsigned int)(index * 1.5);
    unsigned int *new_array = (unsigned int *)realloc(*array, sizeof(unsigned int) * new_length);
    *array = new_array;
    (*array)[index] = value;
    return new_length;
  }

  (*array)[index] = value;

  return length;
}

#ifdef DEBUG

int check_nlist(run_params_t *params, nlist_parameters_t *nlist, double *positions, double rcut)
{
  unsigned n_dims = params->n_dims;
  unsigned n_particles = params->n_particles;
  unsigned int i, j, k, n, ni, nn, n2n, ngn;
  double diff, r;

  for (i = 0, n = 0; i < n_particles; i++)
  {
    //iterate through neighbor list
    nn = 0;
    printf("particle %d nlist: \n", i);
    for (ni = 0; ni < nlist->nlist_count[i]; n++, ni++)
    {
      j = nlist->nlist[n];
      r = 0;

      //distance between particles
      for (k = 0; k < n_dims; k++)
      {
        diff = min_image_dist(positions[j * n_dims + k] - positions[i * n_dims + k], box_size[k]);
        // TODO: remove or not?
        //diff = positions[j * n_dims + k] - positions[i * n_dims + k];
        r += diff * diff;
      }
      if (r < rcut * rcut)
      {
        nn++;
        printf("%d ", j);
      }
    }
    n2n = 0;
    printf("\nparticle %d real neighs: \n", i);
    // stop at i, as we do in nlist code to avoid double couting
    for (j = i + 1; j < n_particles; j++)
    {
      r = 0;
      for (k = 0; k < n_dims; k++)
      {
        diff = positions[j * n_dims + k] - positions[i * n_dims + k];
        r += diff * diff;
      }
      if (r < rcut * rcut)
      {
        n2n++;
        printf("%d ", j);
      }
    }
    ngn = 0;
    printf("\nparticle %d ghost neighs: \n", i);
    for (j = n_particles; j < params->n_ghost_particles + n_particles; j++)
    {
      r = 0;
      for (k = 0; k < n_dims; k++)
      {
        diff = positions[j * n_dims + k] - positions[i * n_dims + k];
        r += diff * diff;
      }
      if (r < rcut * rcut)
      {
        ngn++;
        printf("%d ", j);
      }
    }
    printf("\nParticle %d has %d on nlist and actually has %d real neighbors and %d ghost neighbors\n",
           i, nn, n2n, ngn);
    if (nn != n2n + ngn)
    {
      printf("Invalid nlist detected (see above)\n");
      printf("Nlist: \n");
      n -= nlist->nlist_count[i];
      for (ni = 0; ni < nlist->nlist_count[i]; n++, ni++)
      {
        j = nlist->nlist[n];
        r = 0;

        //distance between particles
        for (k = 0; k < n_dims; k++)
        {
          diff = min_image_dist(positions[j * n_dims + k] - positions[i * n_dims + k], box_size[k]);
          // TODO: remove or not?
          //diff = positions[j * n_dims + k] - positions[i * n_dims + k];
          r += diff * diff;
        }
        printf("%d - %d. r: %f, rcut: %f\n", i, j, sqrt(r), rcut);
      }

      exit(1);
    }
  }
}

#endif