#ifndef ASSIMLSMETHODS_H
#define ASSIMLSMETHODS_H

#include <stdio.h>
#include <stdlib.h>
#include "asynch_interface.h"

extern int np;
extern int my_rank;

typedef struct
{
	unsigned int* fit_states;
	unsigned int* fit_to_universal;
	unsigned int* num_upstream;
	unsigned int** upstream;
	unsigned int dim;
	unsigned int num_fit_states;
} upstream_data;

typedef struct
{
	char* db_filename;
	ConnData* conninfo;	//Wants query to get link ids with gauges, query to download gauge readings
	unsigned int numdata,steps_to_use,least_squares_iters;
	unsigned int* data_locs;
	unsigned int** id_to_assim;
	double inc;
} AssimData;

void ResetSysLS(Link** sys,unsigned int N,UnivVars* GlobalVars,double t_0,double* backup,unsigned int problem_dim,unsigned int num_forcings,TransData* my_data);
double*** ReadSolution(char filename[],unsigned int** id_to_loc,unsigned int N,unsigned int* numlinks,unsigned int** ids,unsigned int** locs,unsigned int** numsteps);
//double*** ReadSolution(char filename[],unsigned int* numlinks,unsigned int** ids,unsigned int** numsteps);
void FindAllDischarges(double*** data,double t,unsigned int numlinks,unsigned int* numsteps,double* d);
unsigned int GaugeDownstream(asynchsolver* asynch,unsigned int** above_gauges,short int** bool_above_gauges,unsigned int* gauges,unsigned int numdata);
int AdjustDischarges_Scale(asynchsolver* asynch,unsigned int* data_locs,double* d,unsigned int numdata,double* x,unsigned int allstates,unsigned int problem_dim);

void Find_Upstream_Links(asynchsolver* asynch,unsigned int problem_dim,short int trim,double inc,unsigned int steps_to_use,unsigned int* data_locs,unsigned int numdata);
void Clean_Upstream_Links(asynchsolver* asynch);
void Free_Upstream_Links(asynchsolver* asynch);

AssimData* Init_AssimData(char* assim_filename,asynchsolver* asynch);
void Free_AssimData(AssimData** assim);
int Download_Gauge_IDs(asynchsolver* asynch,AssimData* Assim);
int GetObservationsDB(AssimData* Assim,unsigned int **id_loc_loc,unsigned int N,unsigned int background_time_unix,double* d);

int ReduceBadDischargeValues(Link** sys,int* assignments,unsigned int N,double* d_full,double* q,unsigned int steps_to_use,unsigned int* data_locs,unsigned int numdata,double* x_start,unsigned int assim_dim,double limit);

int SnapShot_ModelStates(asynchsolver* asynch,unsigned int problem_dim);

int GaugeDataAvailable(AssimData* Assim,unsigned int start_time,unsigned int end_time);

#endif

