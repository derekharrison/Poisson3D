/*
 * export_data.c
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#include <stdio.h>

#include "../inc/user_types.h"


/*-------------------------------------------------------------------------------*/
void export_data(char* file_name,
		         char* type_data,
		         grid_size_t grid_size,
		         double ***ptr)
{
	/*
	 * This function is used to export solver result data to text
	 *
	 * input	grid_size
	 * input	T
	 */

	int i,j,k;
	int nx,ny,nz;

	nx = grid_size.nx;
	ny = grid_size.ny;
	nz = grid_size.nz;

  	FILE *file;
  	file = fopen(file_name,"w");
  	if (file != NULL)
  	{
		fprintf(file, "%s", type_data);

  		for(j=1;j<=ny;j++)
  		{
			fprintf(file,"\n\nj = %i\n\n", j);
			for(k=1;k<=nz;k++)
			{
				for(i=1;i<=nx;i++)
					fprintf(file,"%f\t",ptr[i][j][k]);
				fprintf(file,"\n");
			}
  		}
	 	fclose(file);
  	}
  	else
  	{
		printf("Could not open file");
  	}

}
