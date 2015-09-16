#ifndef _bindingSiteMeasures_h
#define _bindingSiteMeasures_h

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "pdbFile.h"
#include "network.h"
#include "surfaceProbe.h"
#include "modes.h"

void calculateBindingLeverage(pdbFile *macromolecule, modePtr normalModes, float **fnm_t10_data_mode1, float **fnm_t10_data_mode2, float **fnm_t10_data_mode3, float **fnm_t10_data_mode4, float **fnm_t10_data_mode5, float **fnm_t10_data_mode6, float **fnm_t10_data_mode7, float **fnm_t10_data_mode8, float **fnm_t10_data_mode9, float **fnm_t10_data_mode10, int numModes, bindingSite *mergedSites, int numMergedSites, FILE *output, pdbFile *macromolecule_w_heavy);


void map_resid_from_CA_to_heavy(int res_ID, pdbFile *macromolecule_w_heavy); // was not in orig

//float calculateBindingLeverageSite(bindingSite *site, int i, pdbFile *macromolecule, float *R, pdbFile *macromolecule_w_heavy);

float calculateBindingLeverageSite(bindingSite *site, int siteNum, pdbFile *macromolecule, float *R);

float calculateBindingLeverageSite__AddingModes(bindingSite *site, int siteNum, pdbFile *macromolecule, float **fnm_t10_data_modeX);

void calculateLCmergedSite(pdbFile *macromolecule, bindingSite *mergedSites, int numMergedSites, Graph *network, FILE *output);

float calculateLCsite(bindingSite *mergedSites, int i, Graph *network);

void sortArray(float *A, int numMergedSites, int *order);

#endif