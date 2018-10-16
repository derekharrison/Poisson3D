/*
 * poisson_solver_helper.c
 *
 *  Created on: Sep 26, 2018
 *      Author: Derek W. Harrison
 */

#include <math.h>

#include "../inc/memory_functions.h"
#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
void vector_addition(double *v1,
                     double factor_v1,
                     double *v2,
                     double factor_v2,
                     int size_vector,
                     double *output_vector)
{
    /*
     * Addition of vectors of type double
     *
     * input    v1
     * input    factor_v1
     * input    v2
     * input    factor_v2
     * input    size_vector
     * output   output_vector
     */

    int j;

    for (j=1;j<=size_vector;j++)
        output_vector[j] = factor_v1*v1[j] + factor_v2*v2[j];

}


/*-----------------------------------------------------------------------------------------------*/
void dot_product(double *v1,
                 double *v2,
                 int size_vec,
                 double *dot_product)
{
    /*
     * Calculating the dot product of two vectors of type double
     *
     * input    v1
     * input    v2
     * input    size_vec
     * output   dotproduct
     */

    int j;

    double dummy = 0.0;
    for (j = 1; j <= size_vec; j++)
        dummy = dummy + v1[j]*v2[j];

    *dot_product = dummy;

}


/*-----------------------------------------------------------------------------------------------*/
void mat_vec_mult(grid_size_t grid_size,
                  double** A,
                  double* p,
                  double* Ap)
{
    /*
     * Calculating the matrix product A*p where A is the vectorized coefficient matrix
     * and p a vector
     *
     * input    grid_size
     * input    A
     * input    p
     * output   Ap
     */

    int j, nt;

    nt = grid_size.nx*grid_size.ny*grid_size.nz;

    Ap[1] = A[1][4]*p[1] + A[2][3]*p[2] + A[1+grid_size.nx][2]*p[1+grid_size.nx] +
            A[1+grid_size.nx*grid_size.ny][1]*p[1+grid_size.nx*grid_size.ny];

    for (j = 2; j <= grid_size.nx; j++)
        Ap[j] = A[j][3]*p[j-1] + A[j][4]*p[j] + A[j+1][3]*p[j+1] +
        A[j+grid_size.nx][2]*p[j+grid_size.nx] +
        A[j+grid_size.nx*grid_size.ny][1]*p[j+grid_size.nx*grid_size.ny];

    for (j = grid_size.nx + 1; j <= grid_size.nx*grid_size.ny; j++)
        Ap[j] = A[j][2]*p[j-grid_size.nx] + A[j][3]*p[j-1] + A[j][4]*p[j] +
        A[j+1][3]*p[j+1] + A[j+grid_size.nx][2]*p[j+grid_size.nx] +
        A[j+grid_size.nx*grid_size.ny][1]*p[j+grid_size.nx*grid_size.ny];

    for (j = grid_size.nx*grid_size.ny + 1; j <= nt - grid_size.nx*grid_size.ny; j++)
        Ap[j] = A[j][1]*p[j-grid_size.nx*grid_size.ny] + A[j][2]*p[j-grid_size.nx] +
        A[j][3]*p[j-1] + A[j][4]*p[j] + A[j+1][3]*p[j+1] + A[j+grid_size.nx][2]*p[j+grid_size.nx] +
        A[j+grid_size.nx*grid_size.ny][1]*p[j+grid_size.nx*grid_size.ny];

    for (j = nt - grid_size.nx*grid_size.ny + 1; j <= nt - grid_size.nx; j++)
        Ap[j] = A[j][1]*p[j-grid_size.nx*grid_size.ny] + A[j][2]*p[j-grid_size.nx] +
        A[j][3]*p[j-1] + A[j][4]*p[j] + A[j+1][3]*p[j+1] +
        A[j+grid_size.nx][2]*p[j+grid_size.nx];

    for (j = nt - grid_size.nx + 1; j <= nt - 1; j++)
        Ap[j] = A[j][1]*p[j-grid_size.nx*grid_size.ny] + A[j][2]*p[j-grid_size.nx] +
        A[j][3]*p[j-1] + A[j][4]*p[j] + A[j+1][3]*p[j+1];

    Ap[nt] = A[nt][1]*p[nt-grid_size.nx*grid_size.ny] + A[nt][2]*p[nt-grid_size.nx] +
             A[nt][3]*p[nt-1] + A[nt][4]*p[nt];

}


/*-----------------------------------------------------------------------------------------------*/
void generate_coefficient_matrix(domain_size_t domain_size,
                                 fixed_boundaries_t fixed_boundaries,
                                 grid_size_t grid_size,
                                 double gamma,
                                 double (*source)(double x,double y,double z),
                                 double* x_guess,
                                 grid_coordinates_t* grid_coordinates,
                                 double* r,
                                 double** A)
{
    /*
     * Generated the vectorized coefficient matrix for the poisson solver
     *
     * input    domain_size
     * input    grid_size
     * input    fixed_boundaries
     * input    gamma
     * input    source
     * output   grid_coordinates
     * output   r
     * output   A
     */

    double deltax, deltay, deltaz;
    double a1, a2, a3, a4, a5, a6, a7, a8;
    double b1, b2, b3;
    int nn;

    int i,j,k;

    deltax = domain_size.Lx/grid_size.nx;
    deltay = domain_size.Ly/grid_size.ny;
    deltaz = domain_size.Lz/grid_size.nz;

    b1 = -gamma/(deltax*deltax);
    b2 = -gamma/(deltay*deltay);
    b3 = -gamma/(deltaz*deltaz);

    a1 = -3*b1 - 3*b2 - 3*b3;
    a2 = -2*b1 - 3*b2 - 3*b3;
    a3 = -3*b1 - 2*b2 - 3*b3;
    a4 = -2*b1 - 2*b2 - 3*b3;
    a5 = -3*b1 - 3*b2 - 2*b3;
    a6 = -2*b1 - 3*b2 - 2*b3;
    a7 = -3*b1 - 2*b2 - 2*b3;
    a8 = -2*b1 - 2*b2 - 2*b3;

    /*Generating node coordinates*/
    for (i = 1; i <= grid_size.nx; i++)
    for (j = 1; j <= grid_size.ny; j++)
    for (k = 1; k <= grid_size.nz; k++)
    {
        grid_coordinates->X[i][j][k] = i*deltax-deltax/2;
        grid_coordinates->Y[i][j][k] = j*deltay-deltay/2;
        grid_coordinates->Z[i][j][k] = k*deltaz-deltaz/2;
    }

    /*Generating vectorized coefficient matrix*/
    //Generating central coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1]=b3;
        A[nn][2]=b2;
        A[nn][3]=b1;
        A[nn][4]=a8;
        r[nn]= source(grid_coordinates->X[i][j][k],
                      grid_coordinates->Y[i][j][k],
                      grid_coordinates->Z[i][j][k]) -
              (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
               b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] +
               b2*x_guess[nn+grid_size.nx] + b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating corner 1 coefficients
    for (i = 1; i <= 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = 0;
        A[nn][3] = 0;
        A[nn][4] = a1;
        r[nn] = -2*b1*fixed_boundaries.Tw - 2*b2*fixed_boundaries.Ts -
                 2*b3*fixed_boundaries.Tb +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] + b2*x_guess[nn+grid_size.nx] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side a cofficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = a2;
        r[nn] = -2*b2*fixed_boundaries.Ts - 2*b3*fixed_boundaries.Tb +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] +
                 b2*x_guess[nn+grid_size.nx] + b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating corner 2 coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 1; j <= 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = a1;
        r[nn] = -2*b1*fixed_boundaries.Te - 2*b2*fixed_boundaries.Ts -
                 2*b3*fixed_boundaries.Tb +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] + b2*x_guess[nn+grid_size.nx] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side d coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = a3;
        r[nn] = -2*b1*fixed_boundaries.Tw - 2*b3*fixed_boundaries.Tb +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b2*x_guess[nn-grid_size.nx] + A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] +
                 b2*x_guess[nn+grid_size.nx] + b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating face E coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a4;
        r[nn] = -2*b3*fixed_boundaries.Tb +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b2*x_guess[nn-grid_size.nx] + b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] +
                 b1*x_guess[nn+1] + b2*x_guess[nn+grid_size.nx] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side b coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a3;
        r[nn] = -2*b1*fixed_boundaries.Te - 2*b3*fixed_boundaries.Tb +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b2*x_guess[nn-grid_size.nx] + b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] +
                 b2*x_guess[nn+grid_size.nx] + b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating corner 4 coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = a1;
        r[nn] = -2*b1*fixed_boundaries.Tw - 2*b2*fixed_boundaries.Tn -
                 2*b3*fixed_boundaries.Tb +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b2*x_guess[nn-grid_size.nx] + A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side c coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a2;
        r[nn] = -2*b2*fixed_boundaries.Tn - 2*b3*fixed_boundaries.Tb +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b2*x_guess[nn-grid_size.nx] + b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] +
                 b1*x_guess[nn+1] + b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating corner 3 coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a1;
        r[nn] = -2*b1*fixed_boundaries.Te - 2*b2*fixed_boundaries.Tn -
                 2*b3*fixed_boundaries.Tb +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b2*x_guess[nn-grid_size.nx] + b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side e coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = 0;
        A[nn][4] = a5;
        r[nn] = -2*b1*fixed_boundaries.Tw - 2*b2*fixed_boundaries.Ts +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + A[nn][4]*x_guess[nn] +
                 b1*x_guess[nn+1] + b2*x_guess[nn+grid_size.nx] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating face A coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = a6;
        r[nn] = -2*b2*fixed_boundaries.Ts +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b1*x_guess[nn-1] +
                 A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] + b2*x_guess[nn+grid_size.nx] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side f coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 1; j <= 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = a5;
        r[nn] = -2*b1*fixed_boundaries.Te - 2*b2*fixed_boundaries.Ts +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b1*x_guess[nn-1] +
                 A[nn][4]*x_guess[nn] + b2*x_guess[nn+grid_size.nx] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating face D coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = a7;
        r[nn] = -2*b1*fixed_boundaries.Tw +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] + b2*x_guess[nn+grid_size.nx] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating face B coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a7;
        r[nn] = -2*b1*fixed_boundaries.Te +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] + b2*x_guess[nn+grid_size.nx] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side h coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = a5;
        r[nn] = -2*b1*fixed_boundaries.Tw - 2*b2*fixed_boundaries.Tn +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating face C coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a6;
        r[nn] = -2*b2*fixed_boundaries.Tn +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side g coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a5;
        r[nn] = -2*b1*fixed_boundaries.Te - 2*b2*fixed_boundaries.Tn +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] +
                 b3*x_guess[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating corner 5 coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = 0;
        A[nn][4] = a1;
        r[nn] = -2*b1*fixed_boundaries.Tw - 2*b2*fixed_boundaries.Ts -
                 2*b3*fixed_boundaries.Tt +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + A[nn][4]*x_guess[nn] +
                 b1*x_guess[nn+1] + b2*x_guess[nn+grid_size.nx]);
    }

    //Generating side i coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = a2;
        r[nn] = -2*b2*fixed_boundaries.Ts - 2*b3*fixed_boundaries.Tt +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b1*x_guess[nn-1] +
                 A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] + b2*x_guess[nn+grid_size.nx]);
    }

    //Generating corner 6 coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 1; j <= 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = a1;
        r[nn] = -2*b1*fixed_boundaries.Te - 2*b2*fixed_boundaries.Ts -
                 2*b3*fixed_boundaries.Tt +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b1*x_guess[nn-1] +
                 A[nn][4]*x_guess[nn] + b2*x_guess[nn+grid_size.nx]);
    }

    //Generating side l coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = a3;
        r[nn] = -2*b1*fixed_boundaries.Tw - 2*b3*fixed_boundaries.Tt +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] + b2*x_guess[nn+grid_size.nx]);
    }

    //Generating face F coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a4;
        r[nn] = -2*b3*fixed_boundaries.Tt +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] + b1*x_guess[nn+1] +
                 b2*x_guess[nn+grid_size.nx]);
    }

    //Generating side j coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a3;
        r[nn] = -2*b1*fixed_boundaries.Te - 2*b3*fixed_boundaries.Tt +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] + b2*x_guess[nn+grid_size.nx]);
    }

    //Generating corner 8 coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = a1;
        r[nn] = -2*b1*fixed_boundaries.Tw - 2*b2*fixed_boundaries.Tn -
                 2*b3*fixed_boundaries.Tt +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 A[nn][4]*x_guess[nn] + b1*x_guess[nn+1]);
    }

    //Generating side k coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a2;
        r[nn] = -2*b2*fixed_boundaries.Tn - 2*b3*fixed_boundaries.Tt +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 b1*x_guess[nn-1] + A[nn][4]*x_guess[nn] + b1*x_guess[nn+1]);
    }

    //Generating corner 7 coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = a1;
        r[nn] = -2*b1*fixed_boundaries.Te -2*b2*fixed_boundaries.Tn -
                 2*b3*fixed_boundaries.Tt +
                 source(grid_coordinates->X[i][j][k],
                        grid_coordinates->Y[i][j][k],
                        grid_coordinates->Z[i][j][k]) -
                (b3*x_guess[nn-grid_size.nx*grid_size.ny] + b2*x_guess[nn-grid_size.nx] +
                 b1*x_guess[nn-1] + A[nn][4]*x_guess[nn]);
    }

}


/*-----------------------------------------------------------------------------------------------*/
void Ly_solver(grid_size_t grid_size,
               double** L,
               double* r,
               double* y)
{
    /*
     * Solves the linear system Lz = r where L is lower triangular
     *
     * input    grid_size
     * input    L
     * input    r
     * output   y
     */

    int j, nt;

    nt = grid_size.nx*grid_size.ny*grid_size.nz;

    y[1]=r[1]/L[1][4];
    for (j = 2; j <= grid_size.nx; j++)
        y[j]=(r[j]-L[j][3]*y[j-1])/L[j][4];

    for (j = grid_size.nx + 1; j <= grid_size.nx*grid_size.ny; j++)
        y[j]=(r[j] - L[j][2]*y[j-grid_size.nx] - L[j][3]*y[j-1])/L[j][4];

    for (j = grid_size.nx*grid_size.ny + 1; j <= nt; j++)
        y[j]=(r[j] - L[j][1]*y[j-grid_size.nx*grid_size.ny] -
                L[j][2]*y[j-grid_size.nx] - L[j][3]*y[j-1])/L[j][4];

}


/*-----------------------------------------------------------------------------------------------*/
void LTz_solver(grid_size_t grid_size,
                double** L,
                double* y,
                double* z)
{
    /*
     * Solves the linear system L'y = z where L' is upper triangular
     *
     * input    grid_size
     * input    L
     * input    y
     * output   z
     */

    int j, nt;

    nt = grid_size.nx*grid_size.ny*grid_size.nz;

    z[nt] = y[nt]/L[nt][4];
    for(j = nt - 1; j >= nt - grid_size.nx + 1; j--)
        z[j] = (y[j] - L[j+1][3]*z[j+1])/L[j][4];

    for(j = nt - grid_size.nx; j >= nt - grid_size.nx*grid_size.ny + 1; j--)
        z[j] = (y[j] - L[j+1][3]*z[j+1] - L[j+grid_size.nx][2]*z[j+grid_size.nx])/L[j][4];

    for(j = nt - grid_size.nx*grid_size.ny; j >= 1; j--)
        z[j] = (y[j] - L[j+1][3]*z[j+1] - L[j+grid_size.nx][2]*z[j+grid_size.nx] -
                L[j+grid_size.nx*grid_size.ny][1]*z[j+grid_size.nx*grid_size.ny])/L[j][4];

}


/*-----------------------------------------------------------------------------------------------*/
void incomplete_cholesky_factorization(grid_size_t grid_size,
                                       double** A,
                                       double** L)
{
    /*
     * Performs the incomplete cholesky factorization on vectorized matrix A
     *
     * input    grid_size
     * input    A
     * output   L
     */

    int i, j, nt;

    nt = grid_size.nx*grid_size.ny*grid_size.nz;

    L[1][1] = 0;
    L[1][2] = 0;
    L[1][3] = 0;
    L[1][4] = sqrt(A[1][4]);

    for (j = 2; j <= nt; j++)
    for (i = 1; i <= 4; i++)
        if (A[j][i] != 0)
        {
            if (j>=1+grid_size.nx*grid_size.ny)
                L[j][1]=A[j][1]/L[j-grid_size.nx*grid_size.ny][4];
            if (j>=1+grid_size.nx)
                L[j][2]=A[j][2]/L[j-grid_size.nx][4];
            L[j][3]=A[j][3]/L[j-1][4];
            L[j][4]=sqrt(A[j][4] - L[j][1]*L[j][1] - L[j][2]*L[j][2] - L[j][3]*L[j][3]);
        }

}
