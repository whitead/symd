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
 * Naive build without using cells
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
  unsigned int i, j, k;
  double dist, dx;

  for(i = 0; i < n_particles; i++) {
    nlist->nlist_count[i] = 0;
    for(j = i + 1; j < n_particles; j++) {
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
