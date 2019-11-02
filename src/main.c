/*
 * main.c
 *
 *  Created on: September 25 2018
 *      Author: Derek W. Harrison
 *
 *      This code solves the 3D poisson equation, gam*div(grad(T))+q(x,y,z) = 0,
 *      with fixed boundary conditions.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../inc/export_data.h"
#include "../inc/memory_functions.h"
#include "../inc/poisson_solver.h"
#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
static double source_equation(double x,double y,double z)
{
    /*
     * This function specifies the source equation q in the poisson equation
     * gam*div(grad(T))+q(x,y,z) = 0
     *
     * The function can depend on x, y and z coordinates.
     *
     * input    x
     * input    y
     * input    z
     */

    return -sin(M_PI*x) * sin(M_PI*y) * sin(M_PI*z);

}


/*-----------------------------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
    /*
     * The main function of the poisson solver.
     * Parameters and boundary conditions are specified here
     */

    double                         gamma = 0.0;
    bool                      exportData = FALSE;
    domain_size_t            domain_size = {0};
    grid_size_t                grid_size = {0};
    fixed_boundaries_t  fixed_boundaries = {0};
    grid_coordinates_t* grid_coordinates = NULL;
    double***                          T = NULL;

    /*Set parameters and boundary conditions*/
    domain_size.Lx = 1.0;                   //length of domain along x coordinate
    domain_size.Ly = 1.0;                   //length of domain along y coordinate
    domain_size.Lz = 1.0;                   //length of domain along z coordinate

    grid_size.nx = 10;                       //amount of nodes along x coordinate
    grid_size.ny = 10;                       //amount of nodes along y coordinate
    grid_size.nz = 10;                       //amount of nodes along z coordinate

    gamma = 1.0;                            //poisson equation constant

    fixed_boundaries.Tw = 0.0;              //west face boundary condition
    fixed_boundaries.Te = 0.0;              //east face boundary condition
    fixed_boundaries.Ts = 0.0;              //south face boundary condition
    fixed_boundaries.Tn = 0.0;              //north face boundary condition
    fixed_boundaries.Tb = 0.0;              //bottom face boundary condition
    fixed_boundaries.Tt = 0.0;              //top face boundary condition

    exportData = TRUE;                      //data export guard


    /*Allocating memory for output of poisson solver*/
    grid_coordinates = allocate_mem_grid_coordinates(grid_size.nx+1, grid_size.ny+1, grid_size.nz+1);
    T = matrix3D(grid_size.nx+1, grid_size.ny+1, grid_size.nz+1);

    /*Calling poisson solver*/
    poisson_solver(domain_size,
                   fixed_boundaries,
                   grid_size,
                   gamma,
                   source_equation,
                   grid_coordinates,
                   T);

    /*Exporting data*/
    if(exportData)
    {
        export_data("Poisson3D.txt", "Temperature profile", grid_size, T);
        export_data("GridX.txt", "X grid coordinates", grid_size, grid_coordinates->X);
        export_data("GridY.txt", "Y grid coordinates", grid_size, grid_coordinates->Y);
        export_data("GridZ.txt", "Z grid coordinates", grid_size, grid_coordinates->Z);
    }

    /*Deallocate memory*/
    free_grid_coordinates(grid_coordinates, grid_size.nx+1, grid_size.ny+1);
    free_memory_3D(T, grid_size.nx+1, grid_size.ny+1);

    return 0;

}

