/*
 * poisson_solver_helper.h
 *
 *  Created on: Sep 26, 2018
 *      Author: Derek W. Harrison
 */

#ifndef POISSON_SOLVER_HELPER_H_
#define POISSON_SOLVER_HELPER_H_

#include "user_types.h"


void vector_addition(double *v1,
                     double factor_v1,
                     double *v2,
                     double factor_v2,
                     int size_vector,
                     double *output_vector);

void dot_product(double *v1,
                 double *v2,
                 int size_vec,
                 double *dot_product);

void mat_vec_mult(grid_size_t grid_size,
                  double** A,
                  double* p,
                  double* Ap);

void generate_coefficient_matrix(domain_size_t domain_size,
                                 fixed_boundaries_t fixed_boundaries,
                                 grid_size_t grid_size,
                                 double gamma,
                                 double (*source)(double x,double y,double z),
                                 grid_coordinates_t* grid_coordinates,
                                 double* r,
                                 double** A);

void Ly_solver(grid_size_t grid_size,
               double** L,
               double* r,
               double* y);

void LTz_solver(grid_size_t grid_size,
                double** L,
                double* y,
                double* z);

void incomplete_cholesky_factorization(grid_size_t grid_size,
                                       double** A,
                                       double** L);

#endif /* POISSON_SOLVER_HELPER_H_ */
