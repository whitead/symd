#include "nlist.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define COUNT_NLIST

/*
 * Grow and insert a value as needed
 */
unsigned int insert_grow(unsigned int index, unsigned int value, unsigned int** array, unsigned int length);

/*
 * build the neighbor list 
 */
void build_list(double* positions, double* box_size, unsigned int n_dims, unsigned int n_particles, Nlist_Parameters* nlist);

/*
 * Precompute cells used for constructing neighbor lists.
 * Note: this is not ready to be called multiple times. If there
 * ever is barostatting, this needs work.
 */
void build_cells(double* box_size, unsigned int n_dims, unsigned int n_particles, Nlist_Parameters* nlist);


Nlist_Parameters* build_nlist_params(unsigned int n_dims, unsigned int n_particles, double* box_size, double skin, double rcut) {
  
  Nlist_Parameters init = {.rcut = rcut * rcut, .skin = skin * skin , .skin_rcut = (rcut + skin) * (rcut + skin)};
  Nlist_Parameters* nlist = (Nlist_Parameters*) malloc(sizeof(Nlist_Parameters));
  memcpy(nlist, &init, sizeof(Nlist_Parameters));
  nlist->nlist = NULL;
  nlist->last_positions = (double*) malloc(sizeof(double) * n_particles * n_dims);
  build_cells(box_size, n_dims, n_particles, nlist);
  return nlist;
}

void free_nlist(Nlist_Parameters* nlist) {

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


void update_nlist(double* positions, 
	     double* box_size, 
	     unsigned int n_dims, 
	     unsigned int n_particles, 
	     Nlist_Parameters* nlist) {


  if(nlist->nlist == NULL) {
    build_list(positions, box_size, n_dims, n_particles, nlist);
  }
  
  double dist, temp;
  double max1, max2;
  max1 = max2 = 0;
  unsigned int i,k;

  for(i = 0; i < n_particles; i++) {
    dist = 0;
    for(k = 0; k < n_dims; k++) {
      temp = (positions[i * n_dims + k] - nlist->last_positions[i * n_dims + k]);
      temp = min_image_dist(temp, box_size[k]);
      dist += temp * temp;
    }
    if(dist > max1) {
      max1 = dist;
    } else if(dist > max2) {
      max2 = dist;
    }
  }
  if(max1 + max2 > nlist->skin)  {
#ifdef COUNT_NLIST
    printf("updating nlist due to %f + %f > %f\n", max1, max2, nlist->skin);
#endif
    build_list(positions, box_size, n_dims, n_particles, nlist);
  }


  return;
}


unsigned int gen_index_r(int* array, int* dims, unsigned int index, unsigned int cur_dim, unsigned int n_dims, int value, unsigned int output_length) {
  int diff[3] = {0, -1, 1};
  unsigned int i;  
  for(i = 0; index < output_length && i < 3; i++){
    if(cur_dim < n_dims - 1)
      index = gen_index_r(array, dims, index, cur_dim + 1, n_dims, value * dims[cur_dim] + diff[i], output_length);
    else  {
      array[index] = value * dims[cur_dim] + diff[i];
      index++;
    }
  }

  return index;
}


void build_cells(double* box_size, unsigned int n_dims, unsigned int n_particles, Nlist_Parameters* nlist) {

  //rebuild the list
  unsigned int i;
  nlist->cell_number = (int*) malloc(sizeof(int) * n_dims);

  //build cell list
  nlist->cell_number_total = 1;
  double width = sqrt(nlist->skin_rcut);
  for(i = 0; i < n_dims; i++) {
    nlist->cell_number[i] = 
      (int) (box_size[i] + width - 1) / width; // ceiling quotient
    //minimum
    if(nlist->cell_number[i] == 0)
      nlist->cell_number[i] = 1;
    nlist->cell_number[i] += 2; //for ghost cells
    nlist->cell_number_total *= 
      nlist->cell_number[i] - 2; //don't include ghost cells in calculation
  }

#ifdef DEBUG
  printf("cell number: %d\n", nlist->cell_number_total);
#endif //DEBUG

  nlist->ncell_number = (pow(3, n_dims) - 1) / 2 + 1;//number of neighboring cells
  if(nlist->ncell_number > nlist->cell_number_total)//for small systems
    nlist->ncell_number = nlist->cell_number_total;

  nlist->adjacent_cells = (int*) malloc(sizeof(int) * nlist->ncell_number);


  //create offset matrix
  gen_index_r(nlist->adjacent_cells, nlist->cell_number, 0, 0, n_dims, 0, nlist->ncell_number);

#ifdef DEBUG
  printf("offset:\n");
  for(i = 0; i < nlist->ncell_number; i++)
    printf("%d\n", nlist->adjacent_cells[i]);
#endif //DEBUG


  //create mapping for PBCs
  //find how far the cells may extend past real cells with the offset matrix
  nlist->map_offset = 0;
  for(i = 0; i < n_dims; i++)
    nlist->map_offset = nlist->map_offset * nlist->cell_number[i] + 1;
  
  nlist->mapping = 
    (unsigned int*) malloc(sizeof(unsigned int) * (nlist->map_offset * 2 + nlist->cell_number_total));

  //wrapped mapping
  for(i = 0; i < nlist->map_offset; i++) {
    nlist->mapping[i] = 
      nlist->cell_number_total - i % nlist->cell_number_total - 1;
    nlist->mapping[i + nlist->map_offset + nlist->cell_number_total] = 
      i % nlist->cell_number_total;
  }

  //trivial mapping
  for(i = 0; i < nlist->cell_number_total; i++)
    nlist->mapping[i + nlist->map_offset] = i;

#ifdef DEBUG
  printf("mapping:\n");
  for(i = 0; i < nlist->map_offset * 2 + nlist->cell_number_total; i++)
    printf("%d ", nlist->mapping[i]);
  printf("\n");
#endif //DEBUG

  //remove ghost cells from cell count
  for(i = 0; i < n_dims; i++)
    nlist->cell_number[i] -= 2;

  //prepare cell list and head
  nlist->head = (int*) malloc(sizeof(int) * nlist->cell_number_total);
  nlist->cell_list = (int*) malloc(sizeof(int) * n_particles);
  
}

/*
 * Build neighbor list using cells in N-dimensions
 */
void build_list(double* positions, double* box_size, unsigned int n_dims, unsigned int n_particles, Nlist_Parameters* nlist) {

  unsigned int total_n = 0;

  if(nlist->nlist == NULL) {
    //first call, need to set-up some things
    nlist->nlist_count = (unsigned int*) malloc(sizeof(unsigned int) * n_particles);    
    nlist->nlist = (unsigned int*) malloc(sizeof(unsigned int) * n_particles * (n_particles / 2));
    nlist->nlist_length  = n_particles * (n_particles / 2);    
  } 

  //rebuild the list
  unsigned int i, k;
  int j;
  unsigned int ncell;
  double dist, dx;
  int net_dcell;
  int icell;  
  
  //reset cell head
  for(i = 0; i < nlist->cell_number_total; i++)
    nlist->head[i] = -1;


  //fills the list in reverse order. I didn't invent this algorithm. The person who did is a genius.
  for(i = 0; i < n_particles; i++) {
    icell = 0;
    for(k = 0; k < n_dims; k++) {
      icell = (int) (wrap(positions[i * n_dims + k],box_size[k]) / box_size[k]) * nlist->cell_number[k] + icell * nlist->cell_number[k];
    }
    nlist->cell_list[i] = nlist->head[icell];
    nlist->head[icell] = i;      
  }

#ifdef DEBUG
  printf("cell_list:\n");
  for(i = 0; i < n_particles; i++)
    printf("%d ", nlist->cell_list[i]);
  
  printf("\nhead:\n");
  for(i = 0; i < nlist->cell_number_total; i++)
    printf("%d ", nlist->head[i]);
  printf("\n");
#endif //DEBUG

  for(i = 0; i < n_particles; i++) {
    //reset count
    nlist->nlist_count[i] = 0;

    //loop over neighbor cells, with no double counting
    icell = 0;
    //find index of particle using polynomial indexing 
    for(k = 0; k < n_dims; k++)
      icell = (int) (wrap(positions[i * n_dims + k], box_size[k]) / box_size[k]) * nlist->cell_number[k] + icell * nlist->cell_number[k];    
    //loop over neighbor cells
    for(ncell = 0; ncell < nlist->ncell_number; ncell++) {      


      net_dcell = icell + nlist->adjacent_cells[ncell];
      //get head of list, using mapping
      j = nlist->head[nlist->mapping[net_dcell + nlist->map_offset]];

      //-1 marks end of cell list
      while(j != -1) {

	if(nlist->mapping[net_dcell + nlist->map_offset] == 0 && i >= j){
	  //don't double count pairs in self-box
	  break;
	}

	dist = 0;
	for(k = 0; k < n_dims; k++) {
	  dx = positions[i * n_dims + k] - positions[j * n_dims + k];
	  dx = min_image_dist(dx, box_size[k]);
	  dist += dx * dx;
	}

	if(dist < nlist->skin_rcut) {
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
  for(i = 0; i < n_particles; i++) {
    printf("nlist_count[%d]:%d\n", i, nlist->nlist_count[i]);
    for(j = 0; j < nlist->nlist_count[i]; j++)  {
      printf("%d ", nlist->nlist[count]);
      count++;
    }
    printf("\n");
  }
#endif //DEBUG

  //Cache old positions
  for(i = 0; i < n_particles; i++)
    for(k = 0; k < n_dims; k++)
      nlist->last_positions[i * n_dims + k] = positions[i * n_dims + k];
}

/* An insert that grows the array if necessary
 *
 */
unsigned int insert_grow(unsigned int index, unsigned int value, unsigned int** array, unsigned int length) {

  if(index >= length) {
    unsigned int new_length = (unsigned int) (index * 1.5);
    unsigned int* new_array = (unsigned int*) realloc(*array, sizeof(unsigned int) * new_length);
    *array = new_array;
    (*array)[index] = value;
    return new_length;
  }

  (*array)[index] = value;

  return length;
}
