/*
 * export_data.h
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#ifndef EXPORT_DATA_H_
#define EXPORT_DATA_H_

#include "user_types.h"


void export_data(char* file_name,
                 char* type_data,
                 grid_size_t grid_size,
                 double ***ptr);

#endif /* EXPORT_DATA_H_ */
