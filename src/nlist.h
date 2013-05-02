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
  double* last_positions;
  unsigned int* nlist;
  unsigned int* nlist_count;
  unsigned int nlist_length;
  unsigned int* cells;
  unsigned int* cell_i;
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
Nlist_Parameters* build_nlist_params(unsigned int n_dims, unsigned int n_particles, double skin, double rcut);


#endif
