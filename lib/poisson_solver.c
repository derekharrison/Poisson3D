/*
 * poisson_solver.c
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../inc/memory_functions.h"
#include "../inc/poisson_solver_helper.h"
#include "../inc/user_types.h"


/*-------------------------------------------------------------------------------*/
void poisson_solver(domain_size_t domain_size,
					fixed_boundaries_t fixed_boundaries,
					grid_size_t grid_size,
	   			    double gamma,
	   			    double (*source)(double x,double y,double z),
	   			    grid_coordinates_t* grid_coordinates,
	   			    double ***T)
{
	/*
	 * This function solves the 3D poisson equation gamma*div(grad(T))+q(x,y,z,c) = 0
	 * with fixed boundary conditions.
	 *
	 * The equation is solved using a preconditioned conjugate gradient method with
	 * incomplete cholesky factorization as preconditioner.
	 *
	 * Note: a node numbering scheme was used such that the first element of all
	 * arrays start at 1 as opposed to zero. The numbering scheme used in this work is
	 *
	 * nn = i + (j-1)*nx + (k-1)*nx*ny
	 *
	 * where nn is the node number, nx the number of nodes used in the x direction,
	 * ny the number of nodes used in the y direction and i, j, k are indices ranging
	 * from 1 to nt where nt is the total number of nodes in the 3D grid, i.e.: nt = nx*ny*nz
	 * with nz the number of nodes in the z direction.
	 *
	 * input 	domain_size
	 * input 	grid_size
	 * input 	fixed_boundaries
	 * input 	gamma
	 * input 	source(x,y,z)
	 * output	grid_coordinates
	 * output 	T
	 */

	double epsilon, delold, delnew, pAp, error;
	double alpha, B;
	double **A, **L, *Ap, *y, *z, *p, *x, *r;
	int nn, nt, i, j, k, it, imax;

	/*Timing solver*/
	clock_t begin, end;
	double time_spent;
	begin = clock();

	/*Setting variables and coefficients*/
	imax 	= 5000;         //Maximum iterations ICCG
	error 	= 1e-30;        //Tolerance

	nt = grid_size.nx*grid_size.ny*grid_size.nz;

	/*Allocating memory for computation*/
	A 		= matrix2D(nt+1,4+1);
	L 		= matrix2D(nt+1,4+1);
	y 		= matrix1D(nt+1);
	z 		= matrix1D(nt+1);
	p 		= matrix1D(nt+1);
	Ap 		= matrix1D(nt+1);
	x 		= matrix1D(nt+1);
	r 		= matrix1D(nt+1);

	/*Generate coefficient matrix*/
	generate_coefficient_matrix(domain_size,
			    				fixed_boundaries,
			    				grid_size,
			    				gamma,
			    				source,
			    				grid_coordinates,
			    				r,
			    				A);

	/* ICCG preconditioning */
	// Incomplete Cholesky factorization of coefficient matrix A
	incomplete_cholesky_factorization(grid_size, A, L);

	//Solving Ly=r (y = L'z)
	Ly_solver(grid_size, L, r, y);

	//Solving L'z=y
	LTz_solver(grid_size, L, y, z);

	//epsilon = r'*r;
	dot_product(r, r, nt, &epsilon);

	//Setting p = z
	for (j=1;j<=nt;j++)
		p[j]=z[j];

	//Initializing x
	for (j = 1; j <= nt; j++)
		x[j] = 0.0;

    /*Solver iterations*/
	it = 0;
	do
	{
	    //Calculating A*p
		mat_vec_mult(grid_size, A, p, Ap);

	    //delold = r'*z
	    dot_product(r, z, nt, &delold);

	    //p'*Ap
	    dot_product(p, Ap, nt, &pAp);

	    alpha = delold/pAp;

	    //x = x+alpha*p;
	    vector_addition(x, 1.0, p, alpha, nt, x);

	    //r = r-alpha*Ap;
	    vector_addition(r, 1.0, Ap, -alpha, nt, r);

	    //Solving Ly=r (y = L'z)
	    Ly_solver(grid_size, L, r, y);

	    //Solving L'z=y
	    LTz_solver(grid_size, L, y, z);

		//delnew = r'*z
	    dot_product(r, z, nt, &delnew);

	    B = delnew/delold;

	    //p = z + B*p;
    	vector_addition(z, 1.0, p, B, nt, p);

		//r'*r
		dot_product(r, r, nt, &epsilon);

	    epsilon = sqrt(epsilon/nt);
	    it = it + 1;
	}while (it < imax && epsilon > error);

	/*Processing results*/
	for (i=1;i<=grid_size.nx;i++)
	for (j=1;j<=grid_size.ny;j++)
	for (k=1;k<=grid_size.nz;k++)
	{
		nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
		T[i][j][k] = x[nn];
	}

	/*Freeing memory*/
	free_memory_1D(Ap);
	free_memory_1D(y);
	free_memory_1D(z);
	free_memory_1D(p);
	free_memory_1D(x);
	free_memory_1D(r);
	free_memory_2D(A, nt);
	free_memory_2D(L, nt);

	/*Print out some results*/
	end = clock();
	time_spent = (double)(end - begin)/CLOCKS_PER_SEC;

	printf("error: %E\n",epsilon);
	printf("iterations: %d\n", it);
	printf("running time: %f\n", time_spent);
}
