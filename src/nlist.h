#include "min_image.h"

#ifndef NLIST_H_
#define NLIST_H_

/*
 * Implementation of Cell/Verlet neighbor list calculation.
 */

/* When constructing this, the last_positions, skin and rcut must be set.
 * nlist should be set to NULL. skin and rcut are stored as their squares
 *
 */
typedef struct {
  double* last_positions; //used for finding maximum displacement
  unsigned int* nlist; //neighbor list
  unsigned int* nlist_count; // number of neighbors
  int* cell_number; //cells used for constructing neighbor list    
  int* adjacent_cells;//offsets to get adjacent cells
  unsigned int* mapping;//Maps PBC for cells
  int* head;//cell heads. -1 indicates empty
  int* cell_list;//list of cells. -1 indicates terminate.  
  unsigned int nlist_length;    
  unsigned int cell_number_total;//number of cells
  unsigned int ncell_number;//number of adjacent cells
  unsigned int map_offset;//offset for mapping
  const double skin;
  const double rcut;
  const double skin_rcut;
} Nlist_Parameters;


/* Handles knowing when to and actually building/rebuilding neighbor list
 *
 */
void update_nlist(double* positions, double* box_size, unsigned int n_dims, unsigned int n_particles, Nlist_Parameters* nlist);


/* Constructs the structure. skin and rcut should NOT be their sqaures
 *
 */
Nlist_Parameters* build_nlist_params(unsigned int n_dims, unsigned int n_particles, double* box_size, double skin, double rcut);

void free_nlist(Nlist_Parameters* nlist);

#endif
