#ifndef ASSIM_H
#define ASSIM_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <stdlib.h>

#if defined(HAVE_PETSC)

struct Vec;
struct Mat;
struct KSP;

#endif

#include "asynch_interface.h"

//!!!! Mix this with AssimData !!!!
typedef struct
{
    AsynchSolver* asynch;
    double *d_full;
    double *x_start;    //Assimilated initial condition
    double t_b;         //Background time
    double *x_b;        //Background vector
    double *HM_buffer;  //Buffer of the HM Matrix
    Vec rhs;            //Right Hand Side vector
    Vec x;              //Solution of the LS
    Vec B;              //Diagonal of the B Matrix
    Vec R;              //Diagonal of the R Matrix
    Mat HM;
    Mat HTH;
    Mat HMTR;
    KSP ksp;
    Vec invupareas;
    double obs_time_step;
    unsigned int problem_dim;
    unsigned int allstates;
    unsigned int allstates_needed;
    unsigned int num_steps;
    unsigned int *obs_locs;
    unsigned int num_obs;
    unsigned int *above_gauges;
    unsigned int num_above;
    unsigned int assim_dim;
    unsigned int *vareq_shift, *inv_vareq_shift;

    PetscInt *HM_col_indices; //For inserting HM values
    PetscInt *d_indices; //For inserting d values
} Workspace;

typedef struct CustomParams
{
    unsigned int ID, offset, simulation_time_with_data;
} CustomParams;


void Make_Assimilated_Forecasts(AsynchSolver* asynch, unsigned int background_time_unix, double simulation_time_with_data, VEC* backup, Workspace* user, ForecastData* forecaster, AssimData* assim, unsigned int assim_dim, unsigned int forecast_idx);

void MM_mult(double** A, double** B, double** C, unsigned int m, unsigned int inner, unsigned int n);
void VECTOR_Copy(double* u, double* v, unsigned int dim);
int LinearLeastSquares(Workspace* ptr, double* q);
void Print_MATRIX(double** A, unsigned int m, unsigned int n);
void Print_VECTOR(double* v, unsigned int dim);
double*** DownloadGaugeReadings(unsigned int start_time, unsigned int stop_time, unsigned int** id_to_loc, unsigned int N, unsigned int* numlinks, unsigned int** ids, unsigned int** locs, unsigned int** numsteps);
double compute_diff(double* d, double* q, unsigned int size);

int Output_Linkid(double t, VEC y_i, VEC global_params, VEC params, int state, void* user);
int Output_Timestamp(double t, VEC y_i, VEC global_params, VEC params, int state, void* user);

void Init_Output_User_forecastparams(AsynchSolver* asynch);
void Free_Output_User_forecastparams(AsynchSolver* asynch);
void Set_Output_User_forecastparams(AsynchSolver* asynch, unsigned int offset, unsigned int time_with_data);

void Init_Output_PeakflowUser_Offset(AsynchSolver* asynch);
void Free_Output_PeakflowUser_Offset(AsynchSolver* asynch);
void Set_Output_PeakflowUser_Offset(AsynchSolver* asynch, unsigned int offset);


#endif //ASSIM_H
