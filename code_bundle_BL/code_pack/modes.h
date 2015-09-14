#ifndef _modes_h
#define _modes_h

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "pdbFile.h"
#include "network.h"

typedef struct {
	float **fnm_t10_data_mode1;
	float **fnm_t10_data_mode2;
	float **fnm_t10_data_mode3;
	float **fnm_t10_data_mode4;
	float **fnm_t10_data_mode5;
	float **fnm_t10_data_mode6;
	float **fnm_t10_data_mode7;
	float **fnm_t10_data_mode8;
	float **fnm_t10_data_mode9;
	float **fnm_t10_data_mode10;
} modeStr;

typedef modeStr *modePtr;


void read_fnm_t10(float **fnm_t10_data_mode1, float **fnm_t10_data_mode2, float **fnm_t10_data_mode3, float **fnm_t10_data_mode4, float **fnm_t10_data_mode5, float **fnm_t10_data_mode6, float **fnm_t10_data_mode7, float **fnm_t10_data_mode8, float **fnm_t10_data_mode9, float **fnm_t10_data_mode10, FILE *fnm_t10, int n);

#endif