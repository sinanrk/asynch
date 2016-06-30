#ifndef ASSIMMETHODS_H
#define ASSIMMETHODS_H

#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "mathmethods.h"
#include "lapacke.h"
#include "my_scalapack.h"
#include "system.h"
#include "definetype.h"
#include "cblas.h"
#include "clapack.h"

extern int np;
extern int my_rank;

void Diagonalize(MAT* A,MAT* V,VEC* D,int* isuppz);
void Parallel_Diagonalize(MAT* Alocal,MAT* Vlocal,VEC* D,int* desca,int* descz,VEC* work,IVEC* iwork);
void GetState(Link** sys,unsigned int* my_sys,unsigned int N,unsigned int my_N,int* assignments,UnivVars* GlobalVars,double* backup);
void ApplyState(Link** sys,unsigned int N,UnivVars* GlobalVars,double* backup,double backup_time);
double*** ReadSolution(char filename[],unsigned int* numlinks,unsigned int** ids,unsigned int** numsteps);
double FindDischarge(double*** data,unsigned int id,double t,unsigned int numlinks,unsigned int* ids,unsigned int* numsteps,double stddev);
void FindNPRowCol(int* nprow,int* npcol);
double Functional(MAT* B,MAT* R,VEC* d,VEC* analysis,VEC* background);

#endif

