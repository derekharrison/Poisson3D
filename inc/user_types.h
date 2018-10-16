/*
 * user_types.h
 *
 *  Created on: Oct 11, 2018
 *      Author: derek
 */

#ifndef USER_TYPES_H_
#define USER_TYPES_H_


#define TRUE 1
#define FALSE 0

typedef int bool;

typedef struct domain_size_t {
	double Lx;
	double Ly;
	double Lz;
} domain_size_t;

typedef struct fixed_boundaries_t {
	double Tw;
	double Te;
	double Ts;
	double Tn;
	double Tb;
	double Tt;
} fixed_boundaries_t;

typedef struct grid_size_t {
	int nx;
	int ny;
	int nz;
} grid_size_t;

typedef struct grid_coordinates_t {
	double*** X;
	double*** Y;
	double*** Z;
} grid_coordinates_t;

#endif /* USER_TYPES_H_ */
