#ifndef ASSIM_LS_METHODS_H
#define ASSIM_LS_METHODS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <stdlib.h>

#include "asynch_interface.h"
#include "structs.h"

extern int np;
extern int my_rank;

typedef struct UpstreamData
{
    unsigned int* fit_states;       //Holds the index in each state vector of the ith sensitivity at this link.
    unsigned int* fit_to_universal; //Holds universal index of the ith sensitivity at this link.
    unsigned int num_fit_states;    //Number of sensitivity at this link
    unsigned int num_upstreams;     //Number of the upstream links
    Link** upstreams;               //List of the upstream links
    unsigned int num_parents;       //Number of the parents links
    Link** parents;                 //List of the parents links
} UpstreamData;

typedef struct AssimData
{
    char db_filename[ASYNCH_MAX_PATH_LENGTH];
    ConnData conninfo;	        // Query to get link ids with gauges, query to download gauge readings
    unsigned int num_obs;       // Number of observation sites
    unsigned int* obs_locs;     // Link index in the sys[] vector associtated with the site
    unsigned int num_steps;     // Number of time step to use for the optimization
    double obs_time_step;       // Observation time step
    unsigned int max_least_squares_iters;   // Maximum number of LS iterations
    Lookup* id_to_assim;
} AssimData;

void ResetSysLS(Link* sys, unsigned int N, GlobalVars* GlobalVars, double t_0, double* backup, unsigned int problem_dim, unsigned int num_forcings, TransData* my_data);

void FindAllDischarges(double*** data, double t, unsigned int numlinks, unsigned int* numsteps, double* d);
unsigned int GaugeDownstream(const AsynchSolver* asynch, const unsigned int* obs_locs, unsigned int num_obs, unsigned int** above_gauges, bool **is_above_gauges);
int AdjustDischarges(const AsynchSolver* asynch, const unsigned int* obs_locs, const double * obs, unsigned int num_obs, unsigned int problem_dim, double* x);

void FindUpstreamLinks(const AsynchSolver* const asynch, AssimData* const assim, unsigned int problem_dim, bool trim, double obs_time_step, unsigned int num_steps, unsigned int* obs_locs, unsigned int num_obs);
void FindUpstreamLinks2(const AsynchSolver* const asynch, AssimData* const assim, unsigned int problem_dim, bool trim, double obs_time_step, unsigned int num_steps, unsigned int* obs_locs, unsigned int num_obs);

void CleanUpstreamLinks(const AsynchSolver* asynch);
void FreeUpstreamLinks(const AsynchSolver* asynch);

bool InitAssimData(AssimData* assim, const char* assim_filename, AsynchSolver* asynch);
void FreeAssimData(AssimData* assim);

int GetObservationsIds(const AsynchSolver* asynch, AssimData* assim);
int GetObservationsData(const AssimData* assim, const Lookup * const id_loc_loc, unsigned int N, unsigned int background_time_unix, double* d);

bool ReduceBadDischargeValues(Link* sys, int* assignments, unsigned int N, double* d_full, double* q, unsigned int num_steps, unsigned int* data_locs, unsigned int numdata, double* x_start, unsigned int assim_dim, double limit);

int SnapShot_ModelStates(AsynchSolver* asynch, unsigned int problem_dim);

//int GaugeDataAvailable(AssimData* Assim, unsigned int start_time, unsigned int end_time);

#endif //ASSIM_LS_METHODS_H

