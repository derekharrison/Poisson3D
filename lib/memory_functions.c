/*
 * memory_functions.c
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#include <stdlib.h>
#include <stdio.h>

#include "../inc/user_types.h"


/*-------------------------------------------------------------------------------*/
void free_memory_1D_int(int* ptr)
{
	/*
	 * Deallocate 1D array of type int
	 *
	 * input	ptr
	 */

	free(ptr);

}


/*-------------------------------------------------------------------------------*/
void free_memory_1D(double* ptr)
{
	/*
	 * Deallocate 1D array of type double
	 *
	 * input	ptr
	 */

	free(ptr);

}


/*-------------------------------------------------------------------------------*/
void free_memory_2D(double** ptr, int nx)
{
	/*
	 * Deallocate 2D array of type double
	 *
	 * input	ptr
	 * input	nx
	 */

	int i;

	for(i = 0; i < nx; i++)
		free(ptr[i]);

	free(ptr);

}


/*-------------------------------------------------------------------------------*/
void free_memory_3D(double*** ptr, int nx, int ny)
{
	/*
	 * Deallocate 3D array of type double
	 *
	 * input	ptr
	 * input	nx
	 * input	ny
	 */

	int i,j;

	for(i = 0; i < nx; i++)
		for(j = 0; j < ny; j++)
			free(ptr[i][j]);

	for ( i = 0; i < nx; i++)
		free(ptr[i]);

	free(ptr);

}


/*-------------------------------------------------------------------------------*/
int *matrix1D_int(int nx)
{
	/*
	 * Allocate 1D array of type int
	 *
	 * input	nx
	 *
	 * return	m
	 */

    int *m = (int *) calloc(nx, sizeof(int));

    if(!m)
    {
    	printf("Unable to allocate memory for 1D array of type int\n");
    	exit(1);
    }

    return m;

}


/*-------------------------------------------------------------------------------*/
double *matrix1D(int nx)
{
	/*
	 * Allocate 1D array of type double
	 *
	 * input	nx
	 *
	 * return	m
	 */

    double *m = (double *) calloc(nx, sizeof(double));

    if(!m)
    {
    	printf("Unable to allocate memory for 1D array of type double\n");
    	exit(1);
    }

    return m;

}


/*-------------------------------------------------------------------------------*/
double **matrix2D( int nx, int ny)
{
	/*
	 * Allocate 2D array of type double
	 *
	 * input	nx
	 * input	ny
	 *
	 * return	m
	 */

	int i;

	double **m = (double **) calloc ( nx, sizeof( double *));

	if(!m)
	{
		printf("Unable to allocate memory for 2D array of type double\n");
		exit(1);
	}

	for ( i = 0; i < nx; i++)
	{
		m[i] = (double *) calloc ( ny, sizeof( double));

		if(!m[i])
		{
			printf("Unable to allocate memory for 2D array of type double\n");
			exit(2);
		}
	}

	return m;

}


/*-------------------------------------------------------------------------------*/
double ***matrix3D( int nx, int ny, int nz)
{
	/*
	 * Allocate 3D array of type double
	 *
	 * input	nx
	 * input	ny
	 * input	nz
	 *
	 * return	m
	 */

	int i,j;

	double ***m = (double ***) calloc ( nx, sizeof(double **));

	for ( i = 0; i < nx; i++)
	{
		m[i] = (double **) calloc ( ny, sizeof(double *));

		if(!m[i])
		{
			printf("Unable to allocate memory for 3D array of type double\n");
			exit(2);
		}
	}

	for(i = 0; i < nx; i++)
		for(j = 0; j < ny; j++)
		{
			m[i][j] = (double *) calloc(nz, sizeof(double));

			if(!m[i][j])
			{
			 printf("Unable to allocate memory for 3D array of type double\n");
			 exit(3);
			}
		}

	return m;

}


/*-------------------------------------------------------------------------------*/
grid_coordinates_t* allocate_mem_grid_coordinates(int nx, int ny, int nz)
{
	/*
	 * Allocate memory for 3D grid coordinates.
	 *
	 * input	x
	 * input	y
	 * input	z
	 *
	 * return	grid_coordinates
	 */

	grid_coordinates_t* grid_coordinates = malloc(sizeof(grid_coordinates_t));
	grid_coordinates->X					 = matrix3D(nx, ny, nz);
	grid_coordinates->Y					 = matrix3D(nx, ny, nz);
	grid_coordinates->Z					 = matrix3D(nx, ny, nz);

	return grid_coordinates;

}


/*-------------------------------------------------------------------------------*/
void free_grid_coordinates(grid_coordinates_t* grid_coordinates, int nx, int ny)
{
	/*
	 * Deallocate memory used for grid coordinates
	 *
	 * input	grid_coordinates
	 * input 	nx
	 * input	ny
	 */

	free_memory_3D(grid_coordinates->X, nx, ny);
	free_memory_3D(grid_coordinates->Y, nx, ny);
	free_memory_3D(grid_coordinates->Z, nx, ny);
	free(grid_coordinates);

}
