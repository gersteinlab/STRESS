#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "modes.h"


void read_fnm_t10(float **fnm_t10_data_mode1, float **fnm_t10_data_mode2, float **fnm_t10_data_mode3, float **fnm_t10_data_mode4, float **fnm_t10_data_mode5, float **fnm_t10_data_mode6, float **fnm_t10_data_mode7, float **fnm_t10_data_mode8, float **fnm_t10_data_mode9, float **fnm_t10_data_mode10, FILE *fnm_t10, int n)
{
	printf("->read fnm_t10_data \n");
    char line[1024];
	int i;
	int num_atoms = 0;


	while (fgets(line, 1000000, fnm_t10) != NULL)
	{
    	num_atoms++;
	}




	// MODE 1:
	rewind(fnm_t10);
	fgets(line,1024,fnm_t10); 	// gets line "# FourierBasis: cutoff = 2, nbasis = 36"
	fgets(line,1024,fnm_t10);   // gets line "BEGIN MODE *""
	for (i=0; i<n; i++)
	{
		fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode1[i][0]), &(fnm_t10_data_mode1[i][1]), &(fnm_t10_data_mode1[i][2]));
	}
	fgets(line,1024,fnm_t10);  // gets line "END"
	fgets(line,1024,fnm_t10);  // gets line "\n"

	// MODE 2:
	fgets(line,1024,fnm_t10);   // gets line "BEGIN MODE *""
	for (i=0; i<n; i++)
	{
		fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode2[i][0]), &(fnm_t10_data_mode2[i][1]), &(fnm_t10_data_mode2[i][2]));
	}
	fgets(line,1024,fnm_t10);  // gets line "END"
	fgets(line,1024,fnm_t10);  // gets line "\n"

	// MODE 3:
	fgets(line,1024,fnm_t10);   // gets line "BEGIN MODE *""
	for (i=0; i<n; i++)
	{
		fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode3[i][0]), &(fnm_t10_data_mode3[i][1]), &(fnm_t10_data_mode3[i][2]));
	}
	fgets(line,1024,fnm_t10);  // gets line "END"
	fgets(line,1024,fnm_t10);  // gets line "\n"

	// MODE 4:
	fgets(line,1024,fnm_t10);   // gets line "BEGIN MODE *""
	for (i=0; i<n; i++)
	{
		fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode4[i][0]), &(fnm_t10_data_mode4[i][1]), &(fnm_t10_data_mode4[i][2]));
	}
	fgets(line,1024,fnm_t10);  // gets line "END"
	fgets(line,1024,fnm_t10);  // gets line "\n"

	// MODE 5:
	fgets(line,1024,fnm_t10);   // gets line "BEGIN MODE *""
	for (i=0; i<n; i++)
	{
		fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode5[i][0]), &(fnm_t10_data_mode5[i][1]), &(fnm_t10_data_mode5[i][2]));
	}
	fgets(line,1024,fnm_t10);  // gets line "END"
	fgets(line,1024,fnm_t10);  // gets line "\n"

	// MODE 6:
	fgets(line,1024,fnm_t10);   // gets line "BEGIN MODE *""
	for (i=0; i<n; i++)
	{
		fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode6[i][0]), &(fnm_t10_data_mode6[i][1]), &(fnm_t10_data_mode6[i][2]));
	}
	fgets(line,1024,fnm_t10);  // gets line "END"
	fgets(line,1024,fnm_t10);  // gets line "\n"

	// MODE 7:
	fgets(line,1024,fnm_t10);   // gets line "BEGIN MODE *""
	for (i=0; i<n; i++)
	{
		fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode7[i][0]), &(fnm_t10_data_mode7[i][1]), &(fnm_t10_data_mode7[i][2]));
	}
	fgets(line,1024,fnm_t10);  // gets line "END"
	fgets(line,1024,fnm_t10);  // gets line "\n"

	// MODE 8:
	fgets(line,1024,fnm_t10);   // gets line "BEGIN MODE *""
	for (i=0; i<n; i++)
	{
		fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode8[i][0]), &(fnm_t10_data_mode8[i][1]), &(fnm_t10_data_mode8[i][2]));
	}
	fgets(line,1024,fnm_t10);  // gets line "END"
	fgets(line,1024,fnm_t10);  // gets line "\n"

	// MODE 9:
	fgets(line,1024,fnm_t10);   // gets line "BEGIN MODE *""
	for (i=0; i<n; i++)
	{
		fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode9[i][0]), &(fnm_t10_data_mode9[i][1]), &(fnm_t10_data_mode9[i][2]));
	}
	fgets(line,1024,fnm_t10);  // gets line "END"
	fgets(line,1024,fnm_t10);  // gets line "\n"

	// MODE 10:
	fgets(line,1024,fnm_t10);   // gets line "BEGIN MODE *""
	for (i=0; i<n; i++)
	{
		fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode10[i][0]), &(fnm_t10_data_mode10[i][1]), &(fnm_t10_data_mode10[i][2]));
	}
	fgets(line,1024,fnm_t10);  // gets line "END"
	fgets(line,1024,fnm_t10);  // gets line "\n"




	/*
	rewind(fnm_t10);
	fgets(line,1024,fnm_t10); 	//  # FourierBasis: cutoff = 2, nbasis = 36
	for (j=0; j< 10; j++) {
		fgets(line,1024,fnm_t10);  // gets line "BEGIN MODE *""
		//fscanf(fnm_t10, "%1024s", line);
		printf("%s", line);  // prints line "BEGIN MODE *""
		//exit(1);
		for (i=0; i<n; i++)
		{
			fscanf(fnm_t10, "%f %f %f \n", &(fnm_t10_data_mode1[i][0]), &(fnm_t10_data_mode1[i][1]), &(fnm_t10_data_mode1[i][2]));
			printf("%f %f %f \n", fnm_t10_data_mode1[i][0], fnm_t10_data_mode1[i][1], fnm_t10_data_mode1[i][2]);
		}
		fgets(line,1024,fnm_t10);
		printf("%s", line);
		fgets(line,1024,fnm_t10);
		printf("%s", line);
		//exit(1);
	}
	*/




	printf("<-read fnm_t10_data \n");
	return;
}

