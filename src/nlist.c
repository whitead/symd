#include "nlist.h"
#include <string.h>
#include <stdlib.h>


/*
 * Grow and insert a value as needed
 */
void insert_grow(unsigned int index, unsigned int value, unsigned int** array, unsigned int* length);


/*
 * build the neighbor list 
 */
void build_list(double* positions, double* box_size, unsigned int n_dims, unsigned int n_particles, Nlist_Parameters* nlist);


Nlist_Parameters* build_nlist_params(unsigned int n_dims, unsigned int n_particles, double skin, double rcut) {
  
  Nlist_Parameters init = {.rcut = rcut * rcut, .skin = skin * skin , .skin_rcut = (rcut + skin) * (rcut + skin)};
  Nlist_Parameters* nlist = (Nlist_Parameters*) malloc(sizeof(Nlist_Parameters));
  memcpy(nlist, &init, sizeof(Nlist_Parameters));
  nlist->nlist = NULL;
  nlist->last_positions = (double*) malloc(sizeof(double) * n_particles * n_dims);
  return nlist;

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
      dist += temp * temp;
    }
    if(dist > max1) {
      max1 = dist;
    } else if(dist > max2) {
      max2 = dist;
    }
  }
  if(max1 + max2 > nlist->skin)  {
    build_list(positions, box_size, n_dims, n_particles, nlist);
  }


  return;
}


/*
 * Build list using cells
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
  unsigned int cell_number[n_dims];
  unsigned int ncell, cell_number_total, ncell_number;
  double dist, dx;
  int net_dcell;

  ncell_number = 14;
  int neighbors[ncell_number];

  if(n_dims == 3) {
    neighbors[0] = 0 + -1 * cell_number[0] + cell_number[0] * cell_number[1] * 0;
    neighbors[1] = 1 + -1 * cell_number[0] + cell_number[0] * cell_number[1] * 0;
    neighbors[2] = 1 + 0 * cell_number[0] + cell_number[0] * cell_number[1] * 0;
    neighbors[3] = 1 + 1 * cell_number[0] + cell_number[0] * cell_number[1] * 0;
    neighbors[4] = -1 + -1 * cell_number[0] + cell_number[0] * cell_number[1] * -1;
    neighbors[5] = -1 + 0 * cell_number[0] + cell_number[0] * cell_number[1] * -1;
    neighbors[6] = -1 + 1 * cell_number[0] + cell_number[0] * cell_number[1] * -1;
    neighbors[7] = 0 + -1 * cell_number[0] + cell_number[0] * cell_number[1] * -1;
    neighbors[8] = 0 + 0 * cell_number[0] + cell_number[0] * cell_number[1] * -1;
    neighbors[9] = 0 + 1 * cell_number[0] + cell_number[0] * cell_number[1] * -1;
    neighbors[10] = 1 + -1 * cell_number[0] + cell_number[0] * cell_number[1] * -1;
    neighbors[11] = 1 + 0 * cell_number[0] + cell_number[0] * cell_number[1] * -1;
    neighbors[12] = 1 + 1 * cell_number[0] + cell_number[0] * cell_number[1] * -1;
    neighbors[13] = 0 + 0 * cell_number[0] + cell_number[0] * cell_number[1] * 0;
  } else if(n_dims == 2) {
    neighbors[0] = 0 + -1 * cell_number[0];
    neighbors[1] = 1 + -1 * cell_number[0];
    neighbors[2] = 1 + 0 * cell_number[0];
    neighbors[3] = 1 + 1 * cell_number[0];
    neighbors[4] = 0 + 0 * cell_number[0];
    ncell_number = 5;
  } else if(n_dims == 1) {
    neighbors[0] = 0;
    neighbors[0] = 1;
    ncell_number = 2;
  }

  int icell;


  //build cell list
  cell_number_total = 1;
  double width = sqrt(nlist->skin_rcut);
  for(i = 0; i < n_dims; i++) {
    cell_number[i] = (box_size[i] + width - 1) / width; // ceiling quotient
    cell_number_total *= cell_number[i];
  }
  
  int head[cell_number_total];
  int cell_list[n_particles];
  
  for(i = 0; i < cell_number_total; i++)
    head[i] = -1;


  //fills the list in reverse order. I didn't invent this algorithm. The person who did is a genius.
  for(i = 0; i < n_particles; i++) {
    for(k = 0; k < n_dims; k++) {
      icell = (int) ((positions[i * n_dims + k] + box_size[k] / 2) / box_size[k] * cell_number[k]) + icell * cell_number[k];
    }
    cell_list[i] = head[icell];
    head[icell] = i;      
  }


  for(i = 0; i < n_particles; i++) {
    //loop over neighbor cells, with no double counting
    icell = 0;
    //find index of particle using polynomial indexing 
    for(k = 0; k < n_dims; k++)
      icell = (int) ((positions[i * n_dims + k] + box_size[k] / 2) / box_size[k] * cell_number[k]) + icell * cell_number[k];    
    //loop over neighbor cells
    for(ncell = 0; ncell < ncell_number; ncell++) {      

      //get head of list
      net_dcell = icell + neighbors[ncell];
      
      j = head[icell + neighbors[ncell]];
      nlist->nlist_count[i] = 0;

      //-1 marks end of cell list
      while(j != -1) {

	if(j == i)
	  continue;

	dist = 0;
	for(k = 0; k < n_dims; k++) {
	  dx = positions[i * n_dims + k] - positions[j * n_dims + k];
	  dx = min_image_dist(dx, box_size[k]);
	  dist += dx * dx;
	}
	if(dist < nlist->skin_rcut) {
	  insert_grow(total_n, j, &(nlist->nlist), &(nlist->nlist_length));
	  total_n += 1;
	  nlist->nlist_count[i] += 1;
	}

	//get next cell member
	j = cell_list[j];
      }
    }
  }

  //Cache old positions
  for(i = 0; i < n_particles; i++)
    for(k = 0; k < n_dims; k++)
      nlist->last_positions[i * n_dims + k] = positions[i * n_dims + k];
}

/* An insert that grows the array if necessary
 *
 */
void insert_grow(unsigned int index, unsigned int value, unsigned int** array, unsigned int* length) {

  if(index >= *length) {
    unsigned int new_length = (unsigned int) index * 1.5;
    unsigned int* new_array = (unsigned int*) realloc(*array, sizeof(unsigned int) * new_length);
    array = &new_array;
    length = &new_length;
  }

  (*array)[index] = value;

}
