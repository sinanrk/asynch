#ifndef ASSIMLSMETHODS_H
#define ASSIMLSMETHODS_H

#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "mathmethods.h"
#include "system.h"
#include "definetype.h"

extern int np;
extern int my_rank;

void ResetSysLS(Link** sys,unsigned int N,UnivVars* GlobalVars,double t_0,double* backup,unsigned int problem_dim,unsigned int num_forcings,TransData* my_data);
double*** ReadSolution(char filename[],unsigned int** id_to_loc,unsigned int N,unsigned int* numlinks,unsigned int** ids,unsigned int** locs,unsigned int** numsteps);
//double*** ReadSolution(char filename[],unsigned int* numlinks,unsigned int** ids,unsigned int** numsteps);
void FindAllDischarges(double*** data,double t,unsigned int numlinks,unsigned int* numsteps,double* d);
//double FindDischarge(double*** data,unsigned int id,double t,unsigned int numlinks,unsigned int* ids,unsigned int* numsteps);

#endif

