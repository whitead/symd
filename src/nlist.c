#include "nlist.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>



/*
 * Grow and insert a value as needed
 */
unsigned int insert_grow(unsigned int index, unsigned int value, unsigned int** array, unsigned int length);


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
  int cell_number[n_dims];
  unsigned int ncell, cell_number_total, ncell_number;
  double dist, dx;
  int net_dcell;

  //build cell list
  cell_number_total = 1;
  double width = sqrt(nlist->skin_rcut);
  for(i = 0; i < n_dims; i++) {
    cell_number[i] = (int) (box_size[i] + width - 1) / width; // ceiling quotient
    cell_number[i] += 2; //for ghost cells
    cell_number_total *= cell_number[i] - 2; //don't include ghost cells in calculation
  }

#ifdef DEBUG
  printf("cell number: %d\n", cell_number_total);
#endif //DEBUG

  ncell_number = (pow(3, n_dims) - 1) / 2 + 1;//number of neighboring cells
  if(ncell_number > cell_number_total)//for small systems
    ncell_number = cell_number_total;
  int neighbors[ncell_number];


  //create offset matrix
  gen_index_r(neighbors, cell_number, 0, 0, n_dims, 0, ncell_number);

#ifdef DEBUG
  printf("offset:\n");
  for(i = 0; i < ncell_number; i++)
    printf("%d\n", neighbors[i]);
#endif //DEBUG


  //create mapping for PBCs
  //find how far the cells may extend past real cells with the offset matrix
  unsigned int map_offset = 0;
  for(i = 0; i < n_dims; i++)
    map_offset = map_offset * cell_number[i] + 1;
  
  unsigned int mapping[map_offset * 2 + cell_number_total];
  for(i = 0; i < map_offset; i++) {
    mapping[i] = cell_number_total - i % cell_number_total - 1;
    mapping[i + map_offset + cell_number_total] = i % cell_number_total;
  }
  for(i = 0; i < cell_number_total; i++)
    mapping[i + map_offset] = i;

#ifdef DEBUG
  printf("mapping:\n");
  for(i = 0; i < map_offset * 2 + cell_number_total; i++)
    printf("%d ", mapping[i]);
  printf("\n");
#endif //DEBUG

  //remove ghost cells from cell count
  for(i = 0; i < n_dims; i++)
    cell_number[i] -= 2;

  int icell;  
  int head[cell_number_total];
  int cell_list[n_particles];
  
  for(i = 0; i < cell_number_total; i++)
    head[i] = -1;


  //fills the list in reverse order. I didn't invent this algorithm. The person who did is a genius.
  for(i = 0; i < n_particles; i++) {
    icell = 0;
    for(k = 0; k < n_dims; k++) {
      icell = (int) (wrap(positions[i * n_dims + k],box_size[k]) / box_size[k]) * cell_number[k] + icell * cell_number[k];
    }
    cell_list[i] = head[icell];
    head[icell] = i;      
  }

#ifdef DEBUG
  printf("cell_list:\n");
  for(i = 0; i < n_particles; i++)
    printf("%d ", cell_list[i]);
  
  printf("\nhead:\n");
  for(i = 0; i < cell_number_total; i++)
    printf("%d ", head[i]);
  printf("\n");
#endif //DEBUG

  for(i = 0; i < n_particles; i++) {
    //reset count
    nlist->nlist_count[i] = 0;

    //loop over neighbor cells, with no double counting
    icell = 0;
    //find index of particle using polynomial indexing 
    for(k = 0; k < n_dims; k++)
      icell = (int) (wrap(positions[i * n_dims + k], box_size[k]) / box_size[k]) * cell_number[k] + icell * cell_number[k];    
    //loop over neighbor cells
    for(ncell = 0; ncell < ncell_number; ncell++) {      


      net_dcell = icell + neighbors[ncell];
      //get head of list, using mapping
      j = head[mapping[net_dcell + map_offset]];

      //-1 marks end of cell list
      while(j != -1) {

	if(i >= j){
	  //don't double count pairs
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
	j = cell_list[j];
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
