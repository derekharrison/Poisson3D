/*
 * poisson_solver.h
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#ifndef POISSON_SOLVER_H_
#define POISSON_SOLVER_H_

#include "user_types.h"


void poisson_solver(domain_size_t domain_size,
					fixed_boundaries_t fixed_boundaries,
					grid_size_t grid_size,
	   			    double gamma,
	   			    double (*source)(double x,double y,double z),
	   			    grid_coordinates_t* grid_coordinates,
	   			    double ***T);

#endif /* POISSON_SOLVER_H_ */
