/*
 * memory_functions.h
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#ifndef MEMORY_FUNCTIONS_H_
#define MEMORY_FUNCTIONS_H_

#include "user_types.h"


void free_memory_1D_int(int* ptr);
void free_memory_1D(double* ptr);
void free_memory_2D(double** ptr, int nm);
void free_memory_3D(double*** ptr, int nx, int ny);
int *matrix1D_int(int np);
double *matrix1D(int np);
double **matrix2D( int nm, int np);
double ***matrix3D( int nx, int ny, int nz);
grid_coordinates_t* allocate_mem_grid_coordinates(int nx, int ny, int nz);
void free_grid_coordinates(grid_coordinates_t* grid_coordinates, int nx, int ny);

#endif /* MEMORY_FUNCTIONS_H_ */
