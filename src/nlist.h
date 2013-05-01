/*
 * Implementation of Cell/Verlet neighbor list calculation
 */

typedef struct {
  double* last_positions;
  unsigned int* nlist;
  unsigned int* nlist_conut;
  unsigned int* cells;
  unsigned int* cell_i;
} Nlist_Parameters;
