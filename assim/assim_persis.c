#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <libpq-fe.h>
#include <mpi.h>

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#if defined(HAVE_PETSC)
#include <petsc.h>
#endif

#if !defined(_MSC_VER)
#define ASYNCH_SLEEP(seconds) sleep(seconds)
#else
#include <windows.h>
#define ASYNCH_SLEEP(seconds) Sleep((seconds) * 1000)
#endif

#include "assim_ls_methods.h"
#include "asynch_interface.h"
#include "assim_models.h"
#include "forecaster_methods.h"

//!!!! Mix this with AssimData !!!!
typedef struct
{
    double *d_full, *HM_els, t_b, *x_start, *HM_buffer, *x_b;
    Vec *rhs, *x, *invupareas, *B, *R;
    Mat *HM, *HTH, *HMTR;
    AsynchSolver* asynch;
    KSP *ksp;
    double obs_time_step;
    unsigned int problem_dim, allstates, allstates_needed, num_steps;
    unsigned int *obs_locs, num_obs, *above_gauges, num_above, *HM_col_indices, assim_dim, *vareq_shift, *inv_vareq_shift;
    int *d_indices;
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

int Output_Linkid(unsigned int link_id, double t, VEC y_i, VEC global_params, VEC params, int state, void* user);
int Output_Timestamp(unsigned int link_id, double t, VEC y_i, VEC global_params, VEC params, int state, void* user);

void Init_Output_User_forecastparams(AsynchSolver* asynch);
void Free_Output_User_forecastparams(AsynchSolver* asynch);
void Set_Output_User_forecastparams(AsynchSolver* asynch, unsigned int offset, unsigned int time_with_data);

void Init_Output_PeakflowUser_Offset(AsynchSolver* asynch);
void Free_Output_PeakflowUser_Offset(AsynchSolver* asynch);
void Set_Output_PeakflowUser_Offset(AsynchSolver* asynch, unsigned int offset);


int print_out = 1;
int my_rank;
int np;

int main(int argc, char **argv)
{
    //MPI Stuff
    char help[] = "Here is a help message.\n";
    PetscInitialize(&argc, &argv, NULL, help);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //Parse input
    if (argc < 4)
    {
        if (my_rank == 0)
        {
            printf("Command line parameter required:\nA universal variable file (.gbl),\nA forecast file (.fcst),\nA data assimilation file (.das)");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            ASYNCH_SLEEP(10);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    //Declare variables
    PGresult *res;
    char* query = (char*)malloc(1024 * sizeof(char));
    double holder, longest;
    time_t start, stop/*,q_start,q_stop*/;
    unsigned int i, j, k,/*l,m,n,*/current_offset;
    //RKMethod** AllMethods;
    Link* current;
    //char additional[16];	//For output filename
    //srand(time(NULL));	//!!!! Is this needed? !!!!

    AsynchSolver *asynch = Asynch_Init(MPI_COMM_WORLD);

    //Model 15
    //Asynch_Custom_Model(asynch,&SetParamSizes_Assim,&ConvertParams_Assim,&InitRoutines_Assim,&Precalculations_Assim,&ReadInitData_Assim);
    //Model 254
    //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254,&Precalculations_Assim_254,&ReadInitData_Assim_254);
    //Model 254, q
    Asynch_Custom_Model(asynch, &SetParamSizes_Assim_254, &ConvertParams_Assim_254, &InitRoutines_Assim_254_q, &Precalculations_Assim_254, &ReadInitData_Assim_254_q);
    //Model 254, q and s_p
    //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_qsp,&Precalculations_Assim_254,&ReadInitData_Assim_254_qsp);
    //Model 254, q and s_t
    //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_qst,&Precalculations_Assim_254,&ReadInitData_Assim_254_qst);
    //Asynch_Custom_Partitioning(asynch,&Partition_METIS_ByEqs);

    //Data assimilation parameters
/*
    unsigned int num_obs,*ids,*numsteps,*obs_locs;
    double forecast_window = 10.0*24*60;
    double inc = 15.0;
    unsigned int steps_to_use = 10;
    unsigned int least_squares_iters = 1;
*/
//For model 254
//unsigned int problem_dim = 7;	//!!!! Generalize this !!!!
//unsigned int assim_dim = 4;	//!!!! Generalize this !!!!
//For model 254 trim
    unsigned int problem_dim = 4;	//!!!! Generalize this !!!!
    unsigned int assim_dim = 4;	//!!!! Generalize this !!!!

    //Read global file
    if (my_rank == 0)	printf("Reading global file...\n");
    Asynch_Parse_GBL(asynch, argv[1]);
    if (my_rank == 0)	printf("Loading network...\n");
    Asynch_Load_Network(asynch);

    //Read forecast file
    ForecastData* forecaster = Init_ForecastData(argv[2], asynch->globals->string_size);
    if (!forecaster)	MPI_Abort(MPI_COMM_WORLD, 1);
    double forecast_window = forecaster->forecast_window;

    //Read data assimilation file
    AssimData assim;
    InitAssimData(&assim, argv[3], asynch);

    //Make sure the total time is large enough
    holder = Asynch_Get_Total_Simulation_Duration(asynch);
    longest = (holder < forecast_window) ? forecast_window : holder;
    Asynch_Set_Total_Simulation_Duration(asynch, longest);

    //unsigned int obs_locs[] = {2};	//Locations of links with data	!!!! Make sure these are locations and not IDs !!!!
    //unsigned int obs_locs[] = {1,3,4};	//!!!! This should come from ReadSolution. Or should be calculated at least. !!!!
    //obs_locs = (unsigned int*) malloc(3*sizeof(unsigned int));
    //obs_locs[0] = 1; obs_locs[1] = 3; obs_locs[2] = 4;
    //double*** truesolution = ReadSolution("TempData/assim/yobservations.dat",asynch->id_to_loc,N,&num_obs,&ids,&obs_locs,&numsteps);
    //double*** truesolution = ReadSolution("TempData/assim/testobservations.dat",asynch->id_to_loc,asynch->N,&num_obs,&ids,&obs_locs,&numsteps);
    //double*** truesolution = ReadSolution("TempData/assim/254testobservations.dat",asynch->id_to_loc,asynch->N,&num_obs,&ids,&obs_locs,&numsteps);
    //double*** truesolution = DownloadGaugeReadings(1402790400,1405382400,asynch->id_to_loc,N,&num_obs,&ids,&obs_locs,&numsteps);	//June 15 2014 to July 15 2014
    //double*** truesolution = DownloadGaugeReadings(1403654400,1405382400,asynch->id_to_loc,asynch->N,&num_obs,&ids,&obs_locs,&numsteps);	//June 25 2014 to July 15 2014

    //Find the gauged locations
    if (GetObservationsIds(asynch, &assim))
        MPI_Abort(MPI_COMM_WORLD, 1);

    //Finds the link ids upstreams from every gauged locations
    FindUpstreamLinks(asynch, &assim, problem_dim, 1, assim.obs_time_step, assim.num_steps, assim.obs_locs, assim.num_obs);

    if (my_rank == 0)	printf("Partitioning network...\n");
    Asynch_Partition_Network(asynch);
    CleanUpstreamLinks(asynch);
    if (my_rank == 0)	printf("Loading parameters...\n");
    Asynch_Load_Network_Parameters(asynch, 0);
    if (my_rank == 0)	printf("Reading dam and reservoir data...\n");
    Asynch_Load_Dams(asynch);
    if (my_rank == 0)	printf("Setting up numerical error data...\n");
    Asynch_Load_Numerical_Error_Data(asynch);
    if (my_rank == 0)	printf("Initializing model...\n");
    Asynch_Initialize_Model(asynch);
    Setup_Errors(asynch, problem_dim);
    if (my_rank == 0)	printf("Loading initial conditions...\n");
    Asynch_Load_Initial_Conditions(asynch);
    if (my_rank == 0)	printf("Loading forcings...\n");
    Asynch_Load_Forcings(asynch);
    if (my_rank == 0)	printf("Loading output data information...\n");
    Asynch_Load_Save_Lists(asynch);
    if (my_rank == 0)	printf("Finalizing network...\n");
    Asynch_Finalize_Network(asynch);
    if (my_rank == 0)	printf("Calculating initial step sizes...\n");
    Asynch_Calculate_Step_Sizes(asynch);

    Asynch_Set_Total_Simulation_Duration(asynch, holder);

    //Setup output for link id, if needed
    int setup_id = Asynch_Check_Output(asynch, "LinkID");
    int setup_timestamp = Asynch_Check_Output(asynch, "Timestamp");
    if (setup_id || setup_timestamp)
    {
        if (my_rank == 0)	printf("Error: forecaster needs LinkID (%i), Timestamp (%i).\n", setup_id, setup_timestamp);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    Init_Output_User_forecastparams(asynch);
    //Asynch_Set_Output(asynch, "LinkID", ASYNCH_INT, (void(*)(double, VEC, VEC, VEC, int, void*)) &Output_Linkid, NULL, 0);
    Asynch_Set_Output_Int(asynch, "Timestamp", &Output_Timestamp, NULL, 0);
    Set_Output_User_forecastparams(asynch, 0, 0);

    //Prepare output files
    //Asynch_Prepare_Temp_Files(asynch);
    //Asynch_Write_Current_Step(asynch);
    //Asynch_Prepare_Peakflow_Output(asynch);
    //Asynch_Prepare_Output(asynch);
    CreateHaltFile(forecaster->halt_filename);

    //Pull data from asynch
    char dump_filename[ASYNCH_MAX_PATH_LENGTH], filename[ASYNCH_MAX_PATH_LENGTH];
    unsigned int my_N = asynch->my_N, N = asynch->N, *my_sys = asynch->my_sys, **id_to_loc = asynch->id_to_loc, num_forcings = asynch->globals->num_forcings;
    int *assignments = asynch->assignments;
    Link* sys = asynch->sys;
    short int *getting = asynch->getting;
    GlobalVars *globals = asynch->globals;
    Model* custom_model = asynch->custom_model;

    //Set print_time to t_0
    //Asynch_Reset_Temp_Files(asynch,sys[my_sys[0]]->last_t);

    //Find the index of the forcing to use for forecasting
    unsigned int forecast_idx = forecaster->forecasting_forcing;
    if (forecast_idx >= asynch->globals->num_forcings)
    {
        if (my_rank == 0)	printf("Error: No forecasting forcing set.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //Reserve space for backups
    VEC* backup = (VEC*)malloc(N * sizeof(VEC));	//!!!! Same as background x_b? !!!!
    for (i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank || getting[i])
            backup[i] = v_get(sys[i].dim);
        else
            backup[i] = v_get(0);
    }

    //Initialize choices
    unsigned int num_obs = assim.num_obs, *obs_locs = assim.obs_locs, num_steps = assim.num_steps;
    double t_b = 0.0, obs_time_step = assim.obs_time_step;
    double t_f = Asynch_Get_Total_Simulation_Duration(asynch);
    unsigned int dim = assim_dim;
    unsigned int allstates = dim * N;
    //double x_b[allstates];

    // Allocate background
    double  *x_b = malloc(allstates * sizeof(double));
    for (i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank)
        {
            for (j = 0; j < assim_dim; j++)
                x_b[i*assim_dim + j] = sys[i].list->tail->y_approx.ve[j];	//!!!! Need to be able to specify which states are used !!!!
        }
        MPI_Bcast(&(x_b[i*assim_dim]), assim_dim, MPI_DOUBLE, assignments[i], MPI_COMM_WORLD);
    }

    //Other initializations
    //unsigned int numstepstotake;
    unsigned int iterations = (unsigned int)round((t_f - t_b) / obs_time_step);
    double* d = calloc(assim.num_obs, sizeof(double));

    //unsigned int start_idx = 0;
    double* d_full = calloc(num_obs*num_steps, sizeof(double));
    double* x_start = (double*)malloc(allstates * sizeof(double));	//Values used to start asynch solver in tao solvers

    //Find locations unaffected by gauges
    unsigned int *vareq_shift, *inv_vareq_shift;

    //Call model specific data assimilation routines
    //For Model 254
    //unsigned int allstates_needed = Setup_Fitting_Data_Model254(asynch,obs_locs,num_obs);
    //For Model 254 trim, q
    Setup_Fitting_Data_Model254_q(asynch, obs_locs, num_obs);
    //For Model 254 trim, q and s_p
    //Setup_Fitting_Data_Model254_qsp(asynch,obs_locs,num_obs);
    //For Model 254 trim, q and s_t
    //Setup_Fitting_Data_Model254_qst(asynch,obs_locs,num_obs);
    unsigned int allstates_needed = BuildStateShift(asynch, allstates, obs_locs, num_obs, &vareq_shift, &inv_vareq_shift);


    printf("allstates_needed: %u allstates: %u\n", allstates_needed, allstates);

    //For linear least squares
    Vec rhs, x;
    Mat HM, HTH, HMTR;
    KSP ksp;
    double *HM_buffer = (double*)malloc(allstates_needed * sizeof(double));	//!!!! Blah !!!!
    int *HM_col_indices = NULL;
    if (my_rank == 0)
    {
        VecCreateSeq(MPI_COMM_SELF, allstates_needed, &rhs);
        VecCreateSeq(MPI_COMM_SELF, allstates_needed, &x);
        MatCreateSeqDense(MPI_COMM_SELF, num_obs*num_steps, allstates_needed, NULL, &HM);
        MatCreateSeqDense(MPI_COMM_SELF, allstates_needed, allstates_needed, NULL, &HTH);
        MatCreateSeqDense(MPI_COMM_SELF, allstates_needed, num_obs*num_steps, NULL, &HMTR);
        //double* HM_els = (double*) malloc(allstates*sizeof(double));
        HM_col_indices = (int*)malloc(allstates_needed * sizeof(int));
        for (i = 0; i < allstates_needed; i++)	HM_col_indices[i] = i;
        KSPCreate(MPI_COMM_SELF, &ksp);
        KSPSetOperators(ksp, HTH, HTH);
        KSPSetFromOptions(ksp);	//!!!! Do I need this? !!!!
        //MatAssemblyBegin(HTH,MAT_FINAL_ASSEMBLY);
        //MatAssemblyEnd(HTH,MAT_FINAL_ASSEMBLY);
    }

    int* d_indices = (int*)malloc(num_steps*num_obs * sizeof(int));
    for (i = 0; i < num_steps*num_obs; i++)
        d_indices[i] = i;

    //Transfer upstreams areas to all procs
    //char* my_links_needed = (char*) calloc(N,sizeof(char)),*links_needed = (char*) calloc(N,sizeof(char));	//!!!! This seems to cause a problem with Reduce on Helium... !!!!
    short int* my_links_needed = (short int*)calloc(N, sizeof(short int)), *links_needed = (short int*)calloc(N, sizeof(short int));
    double* invupareas = (double*)calloc(N, sizeof(double));
    UpstreamData* updata;
    for (i = 0; i < N; i++)
    {
        //if(assignments[i] == my_rank)	invupareas[i] = 1.0/(sys[i].params.ve[globals->area_idx]*1e3);
        if (assignments[i] == my_rank)	invupareas[i] = 1.0;
        MPI_Bcast(&(invupareas[i]), 1, MPI_DOUBLE, assignments[i], MPI_COMM_WORLD);
    }

    for (i = 0; i < num_obs; i++)	//Links needed for fitting
    {
        if (assignments[obs_locs[i]] == my_rank)
        {
            current = &sys[obs_locs[i]];
            updata = (UpstreamData*)(current->user);
            my_links_needed[obs_locs[i]] = 1;
            //for(j=0;j<current->num_parents;j++)
            //{
            //	for(k=0;k<updata->num_upstreams[j];k++)
            //		my_links_needed[updata->upstreams[j][k]] = 1;
            //}

            for (j = 0; j < updata->num_upstreams; j++)
            {
                my_links_needed[updata->upstreams[j]->location] = 1;
            }
        }
    }

    MPI_Allreduce(my_links_needed, links_needed, N, MPI_SHORT, MPI_LOR, MPI_COMM_WORLD);

    free(my_links_needed);

    //Build weight matrices

    //!!!! Assuming only q is changing !!!!
    unsigned int curr_idx = 0;
    Vec B, R;
    if (my_rank == 0)
    {
        VecCreateSeq(MPI_COMM_SELF, allstates_needed, &B);
        for (i = 0; i < N; i++)
        {
            if (links_needed[i])
                VecSetValue(B, curr_idx++, invupareas[i], INSERT_VALUES);
        }
        VecAssemblyBegin(B);
        VecAssemblyEnd(B);

        VecCreateSeq(MPI_COMM_SELF, num_obs*num_steps, &R);
        for (i = 0; i < num_obs; i++)
        {
            for (j = 0; j < num_steps; j++)
                VecSetValue(R, j*num_obs + i, invupareas[obs_locs[i]] * 10.0, INSERT_VALUES);
        }
        VecAssemblyBegin(R);
        VecAssemblyEnd(R);
    }
    free(links_needed);

    //printf("!!!! Multiplied R by 10. Used 6 m/s (instead of 3). Only used 4 ls iterations. Using tolerance of 10 (instead of 30). Factor is 0.8.!!!!\n");

    /*
        //!!!! Assuming q and s_p are changing !!!!
        unsigned int curr_idx = 0;
        Vec B,R;
        if(my_rank == 0)
        {
            VecCreateSeq(MPI_COMM_SELF,allstates_needed,&B);
            for(i=0;i<N;i++)
            {
                if(links_needed[i])
                {
                    VecSetValue(B,curr_idx++,1.0,INSERT_VALUES);
                    VecSetValue(B,curr_idx++,1e2,INSERT_VALUES);
                }
            }
            VecAssemblyBegin(B);
            VecAssemblyEnd(B);

            VecCreateSeq(MPI_COMM_SELF,num_obs*steps_to_use,&R);
            for(i=0;i<num_obs;i++)
            {
                for(j=0;j<steps_to_use;j++)
                    VecSetValue(R,j*steps_to_use + i,invupareas[obs_locs[i]] * 10.0,INSERT_VALUES);
            }
            VecAssemblyBegin(R);
            VecAssemblyEnd(R);
        }
        free(links_needed);
    */

    //Scale by area
    printf("Scaling init discharges by area...\n");
    //AdjustDischarges_Scale(asynch,obs_locs,d_full,num_obs,x_b,allstates);
    AdjustDischarges(asynch, obs_locs, d_full, num_obs, assim_dim, x_b);

    //Prep PetSC
    if (my_rank == 0)	printf("\nPrepping PetSc...\n");
    Workspace user;
    user.HM = &HM;
    user.HMTR = &HMTR;
    user.B = &B;
    user.R = &R;
    user.HM_els = NULL;
    user.HM_buffer = HM_buffer;
    user.HTH = &HTH;
    user.rhs = &rhs;
    user.HM_col_indices = HM_col_indices;
    user.d_indices = d_indices;
    user.ksp = &ksp;
    //user.d = d;
    user.d_full = d_full;
    user.x_start = x_start;
    user.x = &x;
    user.asynch = asynch;
    user.problem_dim = problem_dim;
    user.assim_dim = assim_dim;
    user.allstates = allstates;
    user.allstates_needed = allstates_needed;
    user.vareq_shift = vareq_shift;
    user.inv_vareq_shift = inv_vareq_shift;
    user.obs_time_step = obs_time_step;
    user.num_steps = assim.num_steps;
    user.obs_locs = obs_locs;
    user.num_obs = num_obs;
    user.t_b = t_b;
    user.x_b = x_b;

    //Print out some information
    unsigned int my_eqs = 0, total_eqs;
    for (i = 0; i < my_N; i++)
        my_eqs += sys[my_sys[i]].dim;
    MPI_Reduce(&my_eqs, &total_eqs, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

    printf("[%i]: Good to go with %u links (%u eqs).\n", my_rank, my_N, my_eqs);
    if (my_rank == 0)
    {
        ASYNCH_SLEEP(1);
        printf("\nNetwork has a total of %u links and %u equations.\n\n", N, total_eqs);
        printf("Making calculations...\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ASYNCH_SLEEP(1);
    start = time(NULL);

    //Begin persistent calculations
    if (my_rank == 0)
        printf("\n\n===================================\nBeginning persistent calculations\n===================================\n");
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    //Make some initializations and checks
    //unsigned int history_time = 5*24*60*60;	//Amount of history to store for hydrographs and peakflows
    short unsigned int hr1 = 0;	//Hour of the day to perform maintainance on database
    short unsigned int hr2 = 12;	//Hour of the day to perform maintainance on database
    unsigned int wait_time = 120;	//Time to ASYNCH_SLEEP if no rainfall data is available
    unsigned int num_tables = 10;
    //unsigned int num_rainsteps = 3;	//Number of rainfall intensities to use for the next forecast
    unsigned int num_rainsteps = forecaster->num_rainsteps;	//Number of rainfall intensities to use for the next forecast
    //if(my_rank == 0 && asynch->globals->increment < num_rainsteps + 3)
    if (my_rank == 0 && asynch->forcings[forecast_idx].increment < num_rainsteps + 3)
        printf("Warning: Increment for rain should probably be %u.\n", num_rainsteps + 3);
    //asynch->forcings[forecast_idx]->increment = num_rainsteps;	//!!!! Not necessary, but makes me feel better. The solvers should really not do the last step where they download nothing. !!!!

    unsigned int /*nextraintime,*/nextforcingtime, background_time_unix;
    short int halt = 0;
    int isnull, repeat_for_errors;
    //int message_buffer[1 + asynch->globals->num_forcings];
    short int vac = 0;	//0 if no vacuum has occured, 1 if vacuum has occured (during a specific hour)
    unsigned int change_time = (unsigned int)asynch->forcings[forecast_idx].file_time * 60 * num_rainsteps;
    unsigned int forecast_time_unix = asynch->forcings[forecast_idx].first_file - change_time;	//Subracting change_time to make the value correct in the loop
    //unsigned int last_file = asynch->forcings[forecast_idx]->last_file;
    //unsigned int first_file = asynch->forcings[forecast_idx]->first_file;
    k = 0;
    for (i = 0; i < N; i++)
        if (backup[i].dim > 0)
            v_copy(asynch->sys[i].list->tail->y_approx, backup[i]);

    //double simulation_time_with_data = asynch->forcings[forecast_idx]->file_time * forecaster->num_rainsteps;
    double simulation_time_with_data = (assim.num_steps - 1) * assim.obs_time_step;
    unsigned int simulation_time_with_data_secs = (int)(simulation_time_with_data + 1e-3) * 60;

    //Setup temp files
    Set_Output_User_forecastparams(asynch, forecast_time_unix, simulation_time_with_data_secs);	//!!!! Should this be done here at all? !!!!
    //Set_Output_PeakflowUser_Offset(asynch,forecast_time_unix);
    Asynch_Set_Total_Simulation_Duration(asynch, forecast_window + simulation_time_with_data);
    Asynch_Prepare_Temp_Files(asynch);

    //Prepare snapshots
    Asynch_Get_Snapshot_Output_Name(asynch, dump_filename);
    dump_filename[strlen(dump_filename) - 4] = '\0';	//Removes .rec	!!!! Uh, is this ok? No chance for memory corruption? !!!!

    //Make some initializations to the database
    if (my_rank == 0 && 1 == 0)
    {
        printf("Making initializations to tables.\n");
        start = time(NULL);

        //Connect to hydrograph database
        ConnectPGDB(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

        //Make sure the hydrographs table exists
        sprintf(query, "SELECT 1 FROM pg_class WHERE relname='%s';", asynch->globals->hydro_table);
        res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT].conn, query);
        if (!PQntuples(res))
        {
            PQclear(res);
            sprintf(query, "CREATE TABLE %s(link_id int,time int,ratio real,discharge real); ALTER TABLE %s SET (autovacuum_enabled = false, toast.autovacuum_enabled = false);", asynch->globals->hydro_table, asynch->globals->hydro_table);
            res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT].conn, query);
            CheckResError(res, "creating hydrographs table");
        }
        else
        {
            PQclear(res);
            sprintf(query, "TRUNCATE %s;", asynch->globals->hydro_table);
            res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT].conn, query);
            CheckResError(res, "truncating hydrographs table");
        }
        PQclear(res);
        DisconnectPGDB(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

        //Make sure the hydrograph tables are set correctly
        CheckPartitionedTable(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT], asynch->globals, forecaster, num_tables, "archive_hydroforecast", "forecast_time", NULL);

        //Clear the future hydrographs in archive
        DeleteFutureValues(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT], num_tables, asynch->globals, "archive_hydroforecast", forecaster->model_name, forecast_time_unix, 1, NULL);
        /*
                sprintf(query,"DELETE FROM hydroforecast_assim WHERE forecast_time >= %u;",forecaster->model_name,forecast_time_unix);
                res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
                CheckResError(res,"deleting future hydroforecasts");
                PQclear(res);

                //Disconnect from hydrograph database, connect to peakflow database
                DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

                ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);

                //Make sure the peakforecast tables are set correctly
                //CheckPeakforecastTable(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT],asynch->globals,forecaster,num_peakflow_tables);

                //Clear all future peakflows
                sprintf(query,"TRUNCATE %s;",asynch->globals->peak_table);
                res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]->conn,query);
                CheckResError(res,"truncating peakforecast table");
                PQclear(res);

                //Disconnect
                DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);
        */
        stop = time(NULL);
        printf("Total time to initialize tables: %.2f.\n", difftime(stop, start));
    }
    else
        printf("!!!! Skipping initialization. Skipping maintenance too. !!!!\n");

    MPI_Barrier(MPI_COMM_WORLD);

    //Start the main loop
    while (!halt)
    {
        if (my_rank == 0)
        {
            time_t now = time(NULL);
            struct tm* now_info = localtime(&now);
            printf("\n\nPass %u\n", k);
            printf("Current time is %s", asctime(now_info));
        }

        //Clear buffers
        Flush_TransData(asynch->my_data);

        //Make some initializations
        //first_file = last_file;
        //last_file = last_file + (unsigned int) asynch->forcings[forecast_idx]->file_time * 60 * num_rainsteps;
        forecast_time_unix += change_time;
        background_time_unix = forecast_time_unix - (assim.num_steps - 1) * (unsigned int)(assim.obs_time_step*60.0 + 1e-3);
        nextforcingtime = forecast_time_unix - 60 * (unsigned int)asynch->forcings[forecast_idx].file_time;	//This is the actual timestamp of the last needed forcing data. This will be downloaded (unlike last_file)

        //Reset each link
        Asynch_Set_System_State(asynch, 0.0, backup);
        Set_Output_User_forecastparams(asynch, forecast_time_unix, simulation_time_with_data_secs);
        //Set_Output_PeakflowUser_Offset(asynch,forecast_time_unix);
        //Asynch_Write_Current_Step(asynch);
        Asynch_Set_Forcing_State(asynch, forecast_idx, 0.0, background_time_unix, forecast_time_unix);

        for (i = 0; i < asynch->globals->num_forcings; i++)	//Set any other database forcings to begin at first_file
        {
            if (asynch->forcings[i].flag == 3)
                Asynch_Set_Forcing_State(asynch, i, 0.0, background_time_unix, asynch->forcings[i].last_file);
        }

        //!!!! Skipped! !!!!
                //Check if a vacuum should be done
                //This will happen at hr1
        //		if(my_rank == 0)	PerformTableMaintainance(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->globals,forecaster,&vac,hr1,num_tables,"archive_hydroforecast");

                //Make sure all buffer flushing is done
        MPI_Barrier(MPI_COMM_WORLD);

        ////Dump data for debugging and recovery
        //if(k % 96 == 0)
        //{
        //	sprintf(filename,"%s%u.rec",dump_filename,forecast_time_unix);
        //	Asynch_Set_Snapshot_Output_Name(asynch,filename);
        //	Asynch_Take_System_Snapshot(asynch,NULL);
        //}

        //Find the next time where rainfall occurs
        do
        {
            if (my_rank == 0)
            {
                ConnectPGDB(forecaster->rainmaps_db);

                //Find the next rainfall time
                time(&start);
                sprintf(query, forecaster->rainmaps_db->queries[0], nextforcingtime);
                res = PQexec(forecaster->rainmaps_db->conn, query);
                CheckResError(res, "checking for new rainfall data");
                time(&stop);
                printf("Total time to check for new rainfall data: %f.\n", difftime(stop, start));
                isnull = PQgetisnull(res, 0, 0);

                PQclear(res);
                DisconnectPGDB(forecaster->rainmaps_db);
            }
            MPI_Bcast(&isnull, 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (isnull)
            {
                if (my_rank == 0)
                {
                    printf("No rainfall values returned from SQL database for forcing %u. %u %u\n", forecast_idx, forecast_time_unix, isnull);
                    PerformTableMaintainance(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT], asynch->globals, forecaster, &vac, hr1, num_tables, "archive_hydroforecast", NULL);
                }

                halt = CheckFinished(forecaster->halt_filename);
                if (halt)
                {
                    sprintf(filename, "%s%u.rec", dump_filename, forecast_time_unix);
                    Asynch_Set_Snapshot_Output_Name(asynch, filename);
                    Asynch_Take_System_Snapshot(asynch, NULL);
                }
                else
                {
                    fflush(stdout);
                    ASYNCH_SLEEP(wait_time);
                }
            }
        } while (isnull && !halt);

        if (halt)	break;

        //Make forecasts
        current_offset = forecast_time_unix;
        Set_Output_User_forecastparams(asynch, current_offset, simulation_time_with_data_secs);	//!!!! Why am I doing this twice? (look like 50 lines up) !!!!
        //Set_Output_PeakflowUser_Offset(asynch,current_offset);
        MPI_Barrier(MPI_COMM_WORLD);
        time(&start);
        Make_Assimilated_Forecasts(asynch, background_time_unix, simulation_time_with_data, backup, &user, forecaster, &assim, assim_dim, forecast_idx);

        MPI_Barrier(MPI_COMM_WORLD);
        if (my_rank == 0)
        {
            time(&stop);
            printf("Time for fitting and forecasts: %.2f\n", difftime(stop, start));
        }
        /*
                //Output some data
                if(my_rank == 0)
                {
                    printf("[%i]: The answer at ID %i at time %.12f is\n",my_rank,asynch->sys[asynch->my_sys[0]]->ID,asynch->sys[asynch->my_sys[0]]->last_t);
                    Print_Vector(asynch->sys[asynch->my_sys[0]]->list->tail->y_approx);
                }
        */
        //Upload the peak data to the database **********************************************************************************************

        //Adjust the table peak
        MPI_Barrier(MPI_COMM_WORLD);
        start = time(NULL);

        repeat_for_errors = Asynch_Create_Peakflows_Output(asynch);
        while (repeat_for_errors > 0)
        {
            if (my_rank == 0)	printf("[%i]: Attempting resend of peakflow data.\n", my_rank);
            ASYNCH_SLEEP(5);
            repeat_for_errors = Asynch_Create_Peakflows_Output(asynch);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if (my_rank == 0)
        {
            stop = time(NULL);
            printf("[%i]: Total time to transfer peak flow data: %.2f\n", my_rank, difftime(stop, start));
        }

        //Upload the hydrographs to the database ********************************************************************************************
        MPI_Barrier(MPI_COMM_WORLD);
        start = time(NULL);

        //Adjust the table hydrographs
        if (my_rank == 0)
        {
            //Make sure database connection is still good
            ConnectPGDB(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

            sprintf(query, "TRUNCATE %s;", asynch->globals->hydro_table);
            res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT].conn, query);
            CheckResError(res, "deleting hydrographs");
            PQclear(res);

            DisconnectPGDB(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        repeat_for_errors = Asynch_Create_Output(asynch, NULL);
        while (repeat_for_errors > 0)
        {
            if (my_rank == 0)	printf("[%i]: Attempting resend of hydrographs data.\n", my_rank);
            ASYNCH_SLEEP(5);
            repeat_for_errors = Asynch_Create_Output(asynch, NULL);
        }

        //Call functions
        if (my_rank == 0)
        {
            //Connect to database
            ConnectPGDB(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

            //Functions for displaying data on IFIS
            if (forecaster->ifis_display)
            {
                //Stage
                repeat_for_errors = 1;
                while (repeat_for_errors)
                {
                    repeat_for_errors = 0;
                    //sprintf(query,"SELECT get_stages_%s();",forecaster->model_name);
                    sprintf(query, "SELECT get_stages_ifc01();");
                    res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT].conn, query);
                    repeat_for_errors = repeat_for_errors || CheckResError(res, "calling stage function");
                    PQclear(res);
                    if (repeat_for_errors)
                    {
                        printf("[%i]: Attempting to call stage function again...\n", my_rank);
                        ASYNCH_SLEEP(5);
                        CheckConnConnection(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
                    }
                }
                /*
                                //Warnings
                                repeat_for_errors = 1;
                                while(repeat_for_errors)
                                {
                                    repeat_for_errors = 0;
                                    sprintf(query,"SELECT update_warnings_%s();",forecaster->model_name);
                                    res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
                                    repeat_for_errors = repeat_for_errors || CheckResError(res,"calling warnings function");
                                    PQclear(res);
                                    if(repeat_for_errors)
                                    {
                                        printf("[%i]: Attempting to call warning function again...\n",my_rank);
                                        ASYNCH_SLEEP(5);
                                        CheckConnConnection(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
                                    }
                                }
                */

                //Run php script
                do
                {
                    repeat_for_errors = system("wget -O /dev/null http://ifisfe.its.uiowa.edu/ifc/php/acquisition/_new_get_ifc_forecast.php");
                    if (repeat_for_errors == -1)
                    {
                        printf("[%i]: Attempting to launch php script again...\n", my_rank);
                        ASYNCH_SLEEP(5);
                    }
                } while (repeat_for_errors == -1);
            }

            //Stage archive
            repeat_for_errors = 0;
            while (repeat_for_errors)
            {
                repeat_for_errors = 0;
                sprintf(query, "ALTER TABLE master_archive_hydroforecast_%s ALTER COLUMN forecast_time SET DEFAULT %u;", forecaster->model_name, current_offset);
                res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT].conn, query);
                repeat_for_errors = repeat_for_errors || CheckResError(res, "setting default value");
                PQclear(res);

                sprintf(query, "SELECT copy_to_archive_hydroforecast_%s();", forecaster->model_name);
                res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT].conn, query);
                repeat_for_errors = repeat_for_errors || CheckResError(res, "calling stage archive function");
                PQclear(res);

                sprintf(query, "ALTER TABLE master_archive_hydroforecast_%s ALTER COLUMN forecast_time DROP DEFAULT;", forecaster->model_name);
                res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT].conn, query);
                repeat_for_errors = repeat_for_errors || CheckResError(res, "dropping default value");
                PQclear(res);

                if (repeat_for_errors)
                {
                    printf("[%i]: Attempting to call stage archive function again...\n", my_rank);
                    ASYNCH_SLEEP(5);
                    CheckConnConnection(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
                }
            }

            //Disconnect
            DisconnectPGDB(&asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
        }

        if (my_rank == 0)
        {
            time(&stop);
            printf("[%i]: Total time to transfer hydrograph data: %.2f\n", my_rank, difftime(stop, start));
        }
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        //Check if program has received a terminate signal **********************************************************************************
        k++;
        halt = CheckFinished(forecaster->halt_filename);

        //If stopping, make a .rec file
        if (halt)
        {
            for (i = 0; i < N; i++)
            {
                current = &asynch->sys[i];
                if (current->list != NULL)
                    v_copy(backup[i], current->list->tail->y_approx);
            }

            sprintf(filename, "%s%u.rec", dump_filename, forecast_time_unix);
            Asynch_Set_Snapshot_Output_Name(asynch, filename);
            Asynch_Take_System_Snapshot(asynch, NULL);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //stop = time(NULL);

    //if(my_rank == 0)	printf("\nTime for calculations: %f. All done!\n",difftime(stop,start));

    //Clean up
    free(query);
    for (i = 0; i < N; i++)
        v_free(&backup[i]);
    free(backup);
    free(d);
    free(d_full);
    free(x_start);
    free(x_start);
    free(vareq_shift);
    free(inv_vareq_shift);

    if (my_rank == 0)
    {
        MatDestroy(&HM);
        MatDestroy(&HTH);
        VecDestroy(&rhs);
        VecDestroy(&x);
        VecDestroy(&R);
        VecDestroy(&B);
        MatDestroy(&HMTR);
    }
    free(HM_buffer);
    free(HM_col_indices);
    free(d_indices);
    if (my_rank == 0)
        KSPDestroy(&ksp);

    free(invupareas);
    FreeAssimData(&assim);
    Free_ForecastData(&forecaster);
    //Free_Output_PeakflowUser_Offset(asynch);
    Free_Output_User_forecastparams(asynch);

    //Petsc clean up
    PetscFinalize();

    //Asynch clean up
    FreeUpstreamLinks(asynch);
    Asynch_Delete_Temporary_Files(asynch);
    Asynch_Free(asynch);

    return 0;
}


void Make_Assimilated_Forecasts(AsynchSolver* asynch, unsigned int background_time_unix, double simulation_time_with_data, VEC* backup, Workspace* user, ForecastData* forecaster, AssimData* assim, unsigned int assim_dim, unsigned int forecast_idx)
{
    time_t q_start, q_stop;
    int *assignments = asynch->assignments;
    short int *getting = asynch->getting;
    unsigned int **id_to_loc = asynch->id_to_loc, N = asynch->N, num_obs = assim->num_obs;
    Link *sys = asynch->sys, *current;
    GlobalVars* globals = asynch->globals;
    unsigned int i, j, l, num_steps = assim->num_steps, assim_window_unix = assim->num_steps * (int)(assim->obs_time_step + 1e-3);
    double t_b = 0.0, obs_time_step = assim->obs_time_step, forecast_window = forecaster->forecast_window, assim_window = assim->num_steps * assim->obs_time_step;
    double *d_full = user->d_full, *x_start = user->x_start, *x_b = user->x_b;
    unsigned int allstates = user->allstates, least_squares_iters = assim->least_squares_iters;
    double *analysis = (double*)malloc(allstates * sizeof(double));	//!!!! Should be removed !!!!
    Model* custom_model = asynch->custom_model;
    //double q[steps_to_use*num_obs];
    double *q = (double*)malloc(num_steps*num_obs * sizeof(double));

    //for(k=0;k<iterations;k++)
    {
        /*
                if(my_rank == 0)
                {
                    printf("\n\n*************************************\n");
                    printf("Iteration %u/%u. Background = %e.\n",k,iterations,t_b);
                }
                MPI_Barrier(MPI_COMM_WORLD);

                steps_to_use = (k < max_steps_to_use) ? k+1 : max_steps_to_use;
        */
        //Set the forecast window
        Asynch_Set_Total_Simulation_Duration(asynch, forecast_window + num_steps*obs_time_step);	//!!!! Is this needed? !!!!

        //Get the new observations
        MPI_Barrier(MPI_COMM_WORLD);
        if (my_rank == 0)
            printf("Downloading observations...\n");
        time(&q_start);

        while (GetObservationsData(assim, id_to_loc, N, background_time_unix, d_full) == -1)
        {
            if (my_rank == 0)	printf("Error downloading observations. Retrying...\n");
            ASYNCH_SLEEP(5);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        /*
                !!!! Need observations !!!!
                GetObseravations(assim,background_time_unix + assim_window_unix,d);	!!!! Only gets newest. Need more for first time through !!!!
                //FindAllDischarges(truesolution,0.0 + (k)*inc,num_obs,numsteps,d);
                for(j=0;j<(steps_to_use-1)*num_obs;j++)
                    d_full[j] = d_full[j+num_obs];
                for(j=0;j<num_obs;j++)
                    d_full[(steps_to_use-1)*num_obs + j] = d[j];
        */
        MPI_Barrier(MPI_COMM_WORLD);
        time(&q_stop);
        if (my_rank == 0)
            printf("Time to get new discharges: %.0f\n", difftime(q_stop, q_start));


        if (print_out && my_rank == 0)
        {
            printf("d_full\n");
            Print_VECTOR(d_full, num_steps*num_obs);
            printf("\n");
        }


        /*
        //if(k == 0)	//!!!! Move this outside the loop. Can I just change minimizer? !!!!
        {
            AdjustDischarges_Scale(asynch,assim->obs_locs,d_full,num_obs,x_b,allstates,assim_dim);
        }
        */


        //Set the initial guess
        for (i = 0; i < N; i++)	//Copy in states that won't be affected by the optimization solver
        {
            for (j = 0; j < assim_dim; j++)
                x_start[i*assim_dim + j] = analysis[i*assim_dim + j] = x_b[i*assim_dim + j];	//!!!! Do I need all these? !!!!
        }


        if (print_out && my_rank == 0)
        {
            //printf("Going in, minimizer...\n");
            //Print_VECTOR(minimizer,allstates_needed);
            printf("Going in, x_start...\n");
            Print_VECTOR(x_start, allstates);

            /*
            double* tempy;
            printf("Guess:\n");
            VecGetArray(Guess,&tempy);
            for(i=0;i<allstates_needed;i++)
            printf("%.15e ",tempy[i]);
            printf("\n");
            VecRestoreArray(Guess,&tempy);
            */
        }

        //Calculate the analysis
        int try_again = 0;
        do
        {
            int iterations = 0;
            double diff = 10.0, error, prev_error = -1.0;
            for (j = 0; j < least_squares_iters; j++)
                //while(diff > 1e-2)
            {
                iterations++;
                LinearLeastSquares(user, q);
                error = compute_diff(d_full, q, num_steps*num_obs);
                if (prev_error >= 0.0)
                {
                    diff = prev_error - error;
                    if (error > prev_error)
                    {
                        if (my_rank == 0)
                        {
                            printf("!!!! LS error got worse. Breaking... !!!!\n");
                            printf("Errors are %f and %f\n", error, prev_error);
                        }

                        //Go back to previous solution
                        for (i = 0; i < allstates; i++)	x_start[i] = analysis[i];

                        break;
                    }
                }
                if (my_rank == 0)	printf("Difference is %f (%f vs %f)\n", diff, error, prev_error);
                prev_error = error;
                for (i = 0; i < allstates; i++)	analysis[i] = x_start[i];
            }
            if (my_rank == 0)
                printf("Total iterations = %i\n", iterations);

            //if(error > 30.0 && !try_again)	//Check if the numerical scheme is having convergence issues
            //if(error > 10.0 && !try_again)
            if (!try_again)
            {
                //try_again = 1;
                try_again = ReduceBadDischargeValues(sys, assignments, N, d_full, q, num_steps, user->obs_locs, num_obs, x_start, assim_dim, 1.0);	//!!!! Not sure what to use for the limit... !!!!
            }
            else	try_again = 0;
        } while (try_again);


        //Copy x_start to analysis  !!!! This shouldn't happen. Try a pointer dance. Or maybe something better... !!!!
        //for(i=0;i<allstates;i++)	analysis[i] = x_start[i];

        if (print_out && my_rank == 0)
        {
            //if(k == 10)
            {
                printf("x_b\n");
                Print_VECTOR(x_b, allstates);
                printf("\n");

                printf("analysis t_b = %e last assim time = %e\n", t_b, t_b + 5.0*(num_steps - 1));
                Print_VECTOR(analysis, allstates);
                //getchar();
            }
        }

        //Switch to model without variational eqs
        //For model 254
        //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Model_254,&Precalculations_Assim_254,&ReadInitData_Assim_254);
        //For model 254 trim
        Asynch_Custom_Model(asynch, &SetParamSizes_Assim_254, &ConvertParams_Assim_254, &InitRoutines_Model_252, &Precalculations_Assim_254, &ReadInitData_Assim_254_q);	//!!!! I think this is ok... !!!!
        Asynch_Initialize_Model(asynch);

        MPI_Barrier(MPI_COMM_WORLD);
        time(&q_start);

        //Advance to get new background
        ResetSysLS(sys, N, globals, t_b, analysis, assim_dim, asynch->globals->num_forcings, asynch->my_data);
        for (i = 0; i < N; i++)				//!!!! Put this into ResetSysLS? !!!!
        {
            if (assignments[i] == my_rank || getting[i])
            {
                current = &sys[i];
                custom_model->initialize_eqs(globals->global_params, current->params, NULL, 0, current->list->head->y_approx, globals->type, current->diff_start, current->no_ini_start, current->user, NULL);
            }
        }

        {
            Asynch_Set_Total_Simulation_Duration(asynch, simulation_time_with_data);
            //globals->maxtime = t_b + inc;
            Asynch_Advance(asynch, 1);

            //Extract background !!!! Isn't there a routine for this? (there needs to be...) !!!!
            for (j = 0; j < N; j++)
            {
                if (assignments[j] == my_rank)
                {
                    for (l = 0; l < assim_dim; l++)
                        x_b[j*assim_dim + l] = backup[j].ve[l] = sys[j].list->tail->y_approx.ve[l];
                }
                MPI_Bcast(&(x_b[j*assim_dim]), assim_dim, MPI_DOUBLE, assignments[j], MPI_COMM_WORLD);
            }
        }


        //Make forecasts
        Asynch_Deactivate_Forcing(asynch, forecast_idx);
        Asynch_Set_Total_Simulation_Duration(asynch, forecast_window + simulation_time_with_data);
        Asynch_Advance(asynch, 1);
        Asynch_Activate_Forcing(asynch, forecast_idx);

        MPI_Barrier(MPI_COMM_WORLD);
        time(&q_stop);
        if (my_rank == 0)	printf("Time for forecast: %.0f\n", difftime(q_stop, q_start));

        //Switch back to model with variational eqs
        //For model 254
        //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254,&Precalculations_Assim_254,&ReadInitData_Assim_254);
        //For model 254 trim, q
        Asynch_Custom_Model(asynch, &SetParamSizes_Assim_254, &ConvertParams_Assim_254, &InitRoutines_Assim_254_q, &Precalculations_Assim_254, &ReadInitData_Assim_254_q);
        //Model 254, q and s_p
        //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_qsp,&Precalculations_Assim_254,&ReadInitData_Assim_254_qsp);
        //Model 254, q and s_t
        //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_qst,&Precalculations_Assim_254,&ReadInitData_Assim_254_qst);
        Asynch_Initialize_Model(asynch);
        /*
                //Go to next time
                t_b += inc;
                user.t_b = t_b;
        */

        /*
                //Reset the temporary files
                if(my_rank == 0)	printf("Creating output file...\n");
        MPI_Barrier(MPI_COMM_WORLD);
        time(&q_start);
                sprintf(additional,"%u",k);	//Add the iteration number to the output files
                Asynch_Create_Output(asynch,additional);
                Asynch_Reset_Temp_Files(asynch,t_b);
        MPI_Barrier(MPI_COMM_WORLD);
        time(&q_stop);
        if(my_rank == 0)	printf("Time to create output: %.0f\n",difftime(q_stop,q_start));
        */
    }

    free(analysis);
    free(q);
}

//This computes the least squares fit assuming the background and analysis difference is linear in the innovations. 
//HM is (num_obs*steps_to_use) X allstates_needed
//HM_els is 1 X allstates (i.e. a 1D array)
//HM_buffer is 1 X allstates_needed
int LinearLeastSquares(Workspace* ptr, double* q)
{
    unsigned int i, j,/*k,m,*/n/*,l,counter*/;
    Link* current;
    //double factor;
    time_t start, stop, start2;
    short int my_link;
    int owner;
    UpstreamData* updata;

    //Unpack ptr
    Workspace* user = (Workspace*)ptr;
    AsynchSolver* asynch = user->asynch;
    unsigned int N = asynch->N;
    Link* sys = asynch->sys;
    GlobalVars* globals = asynch->globals;
    int* assignments = asynch->assignments;
    short int* getting = asynch->getting;
    Mat *HM = user->HM, *HTH = user->HTH, *HMTR = user->HMTR;
    Vec *rhs = user->rhs, d, *x = user->x, *B = user->B, *R = user->R;
    KSP *ksp = user->ksp;
    double *d_els = user->d_full;
    unsigned int *obs_locs = user->obs_locs, assim_dim = user->assim_dim;
    double obs_time_step = user->obs_time_step;
    unsigned int problem_dim = user->problem_dim, allstates = user->allstates, num_steps = user->num_steps, num_obs = user->num_obs;
    int *HM_col_indices = user->HM_col_indices, *d_indices = user->d_indices;
    double t_b = user->t_b, *x_els;
    unsigned int allstates_needed = user->allstates_needed;
    double /**RHS_els,*/*x_start = user->x_start, *HM_buffer = user->HM_buffer;
    unsigned int max_or_steps = num_steps;
    //double q[max_or_steps*num_obs];
    Model* custom_model = asynch->custom_model;
    unsigned int *vareq_shift = user->vareq_shift, *inv_vareq_shift = user->inv_vareq_shift;

    time(&start);

    //Build a vector structure for d !!!! Obviously, this needs to not happen... !!!!
    //if(max_or_steps > allstates_needed)	printf("[%i]: Error: max_or_steps > allstates_needed (%u > %u)\n",my_rank,max_or_steps,allstates_needed);
    VecCreateSeq(MPI_COMM_SELF, max_or_steps*num_obs, &d);	//!!!! This needs to be fixed !!!!
    VecSet(d, 0.0);
    VecSetValues(d, max_or_steps*num_obs, d_indices, d_els, INSERT_VALUES);
    /*
    //!!!! Try to scale the init conditions !!!!
    AdjustDischarges_Scale(asynch,obs_locs,d_els,num_obs,x_start,allstates);
    printf("Adjusted x_start to\n");
    for(i=0;i<allstates;i++)	printf("%.15e ",x_start[i]);
    printf("\n");
    getchar();
    */
    //Initialize the system
    ResetSysLS(sys, N, globals, t_b, x_start, assim_dim, globals->num_forcings, asynch->my_data);
    for (i = 0; i < N; i++)				//!!!! Put this into ResetSysLS? !!!!
        if (assignments[i] == my_rank || getting[i])
            custom_model->initialize_eqs(globals->global_params, sys[i].params, NULL, 0, sys[i].list->head->y_approx, globals->type, sys[i].diff_start, sys[i].no_ini_start, sys[i].user, NULL); //!!!! Should all states be reset? !!!!
            //ReadInitData(globals->global_params,sys[i].params,NULL,0,sys[i].list->head->y_approx,globals->type,sys[i].diff_start,sys[i].no_ini_start,sys[i].user,NULL);	//!!!! Very inefficient. Too many checks. !!!!
    for (i = 0; i < asynch->globals->num_forcings; i++)
    {
        if (asynch->forcings[i].flag == 3)	//!!!! Recurring and binary files need this too !!!!
        {
            //printf("Setting to %u and %u\n",asynch->forcings[i]->first_file,asynch->forcings[i]->last_file);
            Asynch_Set_Forcing_State(asynch, i, t_b, asynch->forcings[i].first_file, asynch->forcings[i].last_file);
        }
    }

    //Advance the system and extract the HM matrix
    for (i = 0; i < max_or_steps; i++)	//!!!! Start at i=1? For i = 0, I don't think we need to set anything... !!!!		//HM here holds the values of M that are needed
    {
        globals->maxtime = t_b + (i)* obs_time_step;
        if (i)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            time(&start2);
            Asynch_Advance(asynch, 0);
            MPI_Barrier(MPI_COMM_WORLD);
            time(&stop);
            if (my_rank == 0)
                printf("Time for advance %u to time %f: %.0f\n", i, globals->maxtime, difftime(stop, start2));
        }


        //printf("ID = %u, t = %e\n",sys[obs_locs[0]]->ID,sys[obs_locs[0]]->last_t);
        //Print_Vector(sys[obs_locs[0]]->list->tail->y_approx);
        //printf("**********\n");

                //Build HM
        for (j = 0; j < num_obs; j++)
        {
            current = &sys[obs_locs[j]];	//!!!! Assumes only discharges !!!!
            owner = assignments[obs_locs[j]];
            my_link = (owner == my_rank);

            //From my link
            if (my_link)
            {
                updata = (UpstreamData*)(current->user);
                //for(n=0;n<allstates;n++)	HM_els[n] = 0.0;

                //Pull out needed data
                for (n = 0; n < allstates_needed; n++)
                    HM_buffer[n] = 0.0;
                for (n = 0; n < updata->num_fit_states; n++)
                {
                    //printf("ID = %u | Loading %e (from %u) into spot %u\n",current->ID,current->list->tail->y_approx.ve[updata->fit_states[n]],updata->fit_states[n],vareq_shift[updata->fit_to_universal[n]]);
                    //ASYNCH_SLEEP(1);
                    HM_buffer[vareq_shift[updata->fit_to_universal[n]]] = current->list->tail->y_approx.ve[updata->fit_states[n]];
                }

                //Extract calculationed q's (Just needed for testing. Maybe...)
                q[i*num_obs + j] = current->list->tail->y_approx.ve[0];
            }

            MPI_Bcast(HM_buffer, allstates_needed, MPI_DOUBLE, owner, MPI_COMM_WORLD);	//!!!! Only proc 0 needs this !!!!
            MPI_Bcast(&(q[i*num_obs + j]), 1, MPI_DOUBLE, owner, MPI_COMM_WORLD);

            unsigned int row_idx = i*num_obs + j;
            //printf("Got %u, %u\n",allstates_needed,updata->num_fit_states);
            //for(n=0;n<allstates_needed;n++)
            //	printf("%e ",HM_buffer[n]);
            //printf("\n");
            //char r = getchar();
            if (my_rank == 0)
                MatSetValues(*HM, 1, &row_idx, allstates_needed, HM_col_indices, HM_buffer, INSERT_VALUES);
        }
    }

    /*
        //Zero out any unused rows
        if(i < steps_to_use)
        {
    printf("Starting at %u, going to %u\n",i*num_obs,steps_to_use*num_obs);
            for(j=0;j<allstates;j++)
                HM_els[j] = 0.0;
            for(i=i*num_obs;i<steps_to_use*num_obs;i++)
                MatSetValues(*HM,1,&i,allstates,cols_allstates_needed,HM_els,INSERT_VALUES);
        }
    */

    //Assemble the HM matrix
    if (my_rank == 0)
    {
        MatAssemblyBegin(*HM, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*HM, MAT_FINAL_ASSEMBLY);

        //MPI_Barrier(MPI_COMM_WORLD);
        time(&start2);
        //Calculate innovations
        double *d_tmp;
        VecGetArray(d, &d_tmp);
        for (i = 0; i < max_or_steps*num_obs; i++)	d_tmp[i] = d_tmp[i] - q[i];
        VecRestoreArray(d, &d_tmp);

        //Build the linear system
        //HMTR is allstates_needed x (num_obs*max_or_steps)
        //HM is (num_obs*max_or_steps) x allstates_needed
        MatTranspose(*HM, MAT_REUSE_MATRIX, HMTR);
        MatDiagonalScale(*HMTR, NULL, *R);
        MatMatMult(*HMTR, *HM, MAT_REUSE_MATRIX, PETSC_DEFAULT, HTH);
        MatDiagonalSet(*HTH, *B, ADD_VALUES);
        /*
                MatTransposeMatMult(*HM,*HM,MAT_REUSE_MATRIX,PETSC_DEFAULT,HTH);
                for(i=0;i<allstates_needed;i++)
                {
                    //MatSetValue(*HTH,i,i,1.0,ADD_VALUES);  //!!!! To skip hillslope !!!!
                    //if(i%2)	MatSetValue(*HTH,i,i,1e3,ADD_VALUES);	//Used for s_p, I think...
                    //if(i%2)	MatSetValue(*HTH,i,i,1e2,ADD_VALUES);	//Used for s_t
                    //else	MatSetValue(*HTH,i,i,1.0,ADD_VALUES);
                }
        */
        MatMult(*HMTR, d, *rhs);
        MatAssemblyBegin(*HTH, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*HTH, MAT_FINAL_ASSEMBLY);

        //MPI_Barrier(MPI_COMM_WORLD);
        time(&stop);
        if (my_rank == 0)
            printf("Time for matrix computations: %.0f\n", difftime(stop, start2));

        //Compute analysis
    //MPI_Barrier(MPI_COMM_WORLD);
        time(&start2);
        KSPSolve(*ksp, *rhs, *x);
        //MPI_Barrier(MPI_COMM_WORLD);
        time(&stop);
        if (my_rank == 0)
            printf("Time for inversion: %.0f\n", difftime(stop, start2));

        //Copy new solution to x_start
        VecGetArray(*x, &x_els);
        //for(i=0;i<num_above;i++)
        for (i = 0; i < allstates_needed; i++)	//!!!! I think this is right... !!!!
        {
            //printf("i = %u %u\n",i,inv_vareq_shift[i]);
            //ASYNCH_SLEEP(1);
            x_start[inv_vareq_shift[i]] += x_els[i];
            //x_start[above_gauges[i]*assim_dim] += x_els[i];	//!!!! To skip hillslope !!!!
            //for(j=0;j<assim_dim;j++)
            //	x_start[above_gauges[i]*assim_dim+j] += x_els[i*assim_dim+j];
        }
        VecRestoreArray(*x, &x_els);
    }

    //Send solution to everyone
    MPI_Bcast(x_start, allstates, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //if(print_out && my_rank == 0)
    //{
    //unsigned int idxm[num_obs*steps_to_use];
    //double temp_matptr[(num_obs*steps_to_use*allstates_needed > allstates_needed*allstates_needed) ? num_obs*steps_to_use*allstates_needed : allstates_needed*allstates_needed];
    //for(i=0;i<num_obs*steps_to_use;i++)	idxm[i] = i;
    //
    //printf("x_start\n");
    //for(i=0;i<allstates;i++)	printf("%.15e ",x_start[i]);
    //printf("\n");
    //
    //double* temp_ptr;
    //printf("difference (x)\n");
    //VecGetArray(*x,&temp_ptr);
    //for(i=0;i<allstates_needed;i++)	printf("%.15e ",temp_ptr[i]);
    //printf("\n");
    //VecRestoreArray(*x,&temp_ptr);
    //
    //printf("d\n");
    //VecGetArray(d,&temp_ptr);
    //for(i=0;i<num_obs*max_or_steps;i++)	printf("%.15e ",temp_ptr[i]);
    //printf("\n");
    //
    //printf("HM\n");
    //MatGetValues(*HM,num_obs*max_or_steps,idxm,allstates_needed,cols_allstates_needed,temp_matptr);
    //for(i=0;i<num_obs*max_or_steps;i++)
    //{
    //	for(j=0;j<allstates_needed;j++)
    //		printf("%.15e ",temp_matptr[i*allstates_needed + j]);
    //	printf(";\n");
    //}
    //
    //printf("HTH\n");
    //MatGetValues(*HTH,allstates_needed,cols_allstates_needed,allstates_needed,cols_allstates_needed,temp_matptr);
    //for(i=0;i<allstates_needed;i++)
    //{
    //	for(j=0;j<allstates_needed;j++)
    //		printf("%.15e ",temp_matptr[i*allstates_needed + j]);
    //	printf(";\n");
    //}
    //
    //}

    if (print_out)
    {
        //Get q's produced from analysis (for testing)
        if (my_rank == 0)
        {
            //printf("q before\n");
            //Print_VECTOR(q,max_or_steps*num_obs);
        }
        double first_diff = compute_diff(d_els, q, max_or_steps*num_obs);

        ResetSysLS(sys, N, globals, t_b, x_start, problem_dim, globals->num_forcings, asynch->my_data);
        for (i = 0; i < N; i++)
            if (assignments[i] == my_rank || getting[i])
                //ReadInitData(globals->global_params,sys[i].params,NULL,0,sys[i].list->head->y_approx,globals->type,sys[i].diff_start,sys[i].no_ini_start,sys[i].user,NULL);
                custom_model->initialize_eqs(globals->global_params, sys[i].params, NULL, 0, sys[i].list->head->y_approx, globals->type, sys[i].diff_start, sys[i].no_ini_start, sys[i].user, NULL); //!!!! Should all states be reset? !!!!
        for (i = 0; i < asynch->globals->num_forcings; i++)
        {
            if (asynch->forcings[i].flag == 3)	//!!!! I think .mon and binary files need this too !!!!
                Asynch_Set_Forcing_State(asynch, i, t_b, asynch->forcings[i].first_file, asynch->forcings[i].last_file);
        }
        for (i = 0; i < max_or_steps; i++)
        {
            globals->maxtime = t_b + (i)* obs_time_step;
            if (i)	Asynch_Advance(asynch, 0);

            for (j = 0; j < num_obs; j++)
            {
                owner = assignments[obs_locs[j]];
                my_link = (owner == my_rank);

                if (my_link)
                    q[i*num_obs + j] = sys[obs_locs[j]].list->tail->y_approx.ve[0];

                MPI_Bcast(&(q[i*num_obs + j]), 1, MPI_DOUBLE, owner, MPI_COMM_WORLD);
            }
        }

        if (my_rank == 0)
        {
            //printf("q after\n");
            //Print_VECTOR(q,max_or_steps*num_obs);

            double second_diff = compute_diff(d_els, q, max_or_steps*num_obs);
            printf("\nDifferences between q and data are %e %e\n\n", first_diff, second_diff);
        }
    }



    //Clean up
    VecDestroy(&d);	//!!!! Blah !!!!

    MPI_Barrier(MPI_COMM_WORLD);
    time(&stop);
    if (my_rank == 0)	printf("Total time for linear least squares fit: %.0f\n", difftime(stop, start));

    return 0;
    //if(second_diff < first_diff)	return 1;
    //else				return 0;
}


double compute_diff(double* d, double* q, unsigned int size)
{
    unsigned int i;
    double result = 0.0;

    for (i = 0; i < size; i++)
        result += (d[i] - q[i]) * (d[i] - q[i]);

    return pow(result, 0.5);
}


//v = u
void VECTOR_Copy(double* u, double* v, unsigned int dim)
{
    unsigned int i;
    for (i = 0; i < dim; i++)
        v[i] = u[i];
}

//C = A*B
//A= m x inner, B= inner x p, C= m x n
void MM_mult(double** A, double** B, double** C, unsigned int m, unsigned int inner, unsigned int n)
{
    unsigned int i, j, k;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)	C[i][j] = 0.0;
        for (k = 0; k < inner; k++)
        {
            for (j = 0; j < n; j++)
                C[i][j] += A[i][k] * B[k][j];
        }
    }
}

void Print_MATRIX(double** A, unsigned int m, unsigned int n)
{
    unsigned int i, j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)	printf("%.15e ", A[i][j]);
        printf(";\n");
    }
    printf("\n");
}

void Print_VECTOR(double* v, unsigned int dim)
{
    unsigned int i;

    for (i = 0; i < dim; i++)
        printf("%.15e ", v[i]);
    printf(";\n");
}


//Read into memory the times and discharges stored in a .dat file.
double*** DownloadGaugeReadings(unsigned int start_time, unsigned int stop_time, unsigned int** id_to_loc, unsigned int N, unsigned int* numlinks, unsigned int** ids, unsigned int** locs, unsigned int** numsteps)
{
    unsigned int i, j/*,k*/;
    double ***data = NULL;
    char query[1028];

    if (my_rank == 0)
    {
        printf("Downloading gauge data...\n");
        PGresult *res;

        //Hard coding!
        //unsigned int outlet = 434478;  //Turkey River above French Hollow
        unsigned int outlet = 434514;	//Turkey River at Garber
        //unsigned int outlet = 307864;  //Half Squaw Creek
        //unsigned int outlet = 292254;	//Squaw Creek at Ames
        //ConnData* conninfo = CreateConnData("dbname=model_ifc host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
        ConnData conninfo;
        ConnData_Init(&conninfo, "dbname=arch_usgs host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
        ConnectPGDB(&conninfo);

        //Get link ids of gauges
        //sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
		//  SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, sensor_data AS C \
		//  WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND B.link_id = A.link_id AND B.link_id != 311903 ORDER BY A.link_id;",outlet);
        //sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
		//	SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, sensor_data AS C \
		//	WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND B.link_id = A.link_id AND B.link_id != 311903 AND B.link_id != 301218 AND B.link_id != 305680 ORDER BY A.link_id;",outlet);
        sprintf(query, "WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
			SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, master_usgs_gauges AS C \
			WHERE B.link_id = A.link_id AND A.id = C.ifis_id AND A.link_id != 421097 AND A.link_id != 434582 ORDER BY link_id;", outlet);
        res = PQexec(conninfo.conn, query);
        CheckResError(res, "getting list of bridge sensor link ids");
        *numlinks = PQntuples(res);

        //Allocate space
        *ids = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
        *locs = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
        *numsteps = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
        data = (double***)malloc(*numlinks * sizeof(double**));

        //Store the ids
        for (i = 0; i < *numlinks; i++)
            (*ids)[i] = atoi(PQgetvalue(res, i, 0));

        //Download data
        for (i = 0; i < *numlinks; i++)
        {
            //sprintf(query,"SELECT (extract('epoch' FROM date_trunc('minute', C.time2)) - %u)/60 AS unix_time, getdischarges((dist_bottom - measured_dist)*0.0328084::real,A.id)*0.0283168 AS q \
			//	FROM env_pois_adv AS A, sensor_data AS C \
			//	WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND A.link_id = %u AND to_timestamp(%u) <= C.time2 AND C.time2 <= to_timestamp(%u) \
			//	ORDER BY unix_time;",start_time,(*ids)[i],start_time,stop_time);
            sprintf(query, "SELECT (unix_timestamp-%u)/60 AS t,discharge*0.0283168 AS q FROM env_pois_adv AS A, master_usgs_gauges AS B \
				WHERE A.id = B.ifis_id AND A.link_id = %u AND %u <= B.unix_timestamp AND B.unix_timestamp <= %u ORDER BY t;", start_time, (*ids)[i], start_time, stop_time);
            res = PQexec(conninfo.conn, query);
            CheckResError(res, "downloading bridge sensor data");
            (*numsteps)[i] = PQntuples(res);
            data[i] = (double**)malloc((*numsteps)[i] * sizeof(double*));

            for (j = 0; j < (*numsteps)[i]; j++)
            {
                data[i][j] = (double*)malloc(2 * sizeof(double));
                data[i][j][0] = atof(PQgetvalue(res, j, 0));
                data[i][j][1] = atof(PQgetvalue(res, j, 1));
            }
        }

        //Find locations from ids
        for (i = 0; i < *numlinks; i++)
            (*locs)[i] = find_link_by_idtoloc((*ids)[i], id_to_loc, N);

        //Cleanup
        PQclear(res);
        DisconnectPGDB(&conninfo);
        ConnData_Free(&conninfo);
    }

    /*
        //Send data to all procs
        //!!!! Would be better if only data for links assigned to this proc were available !!!!
        //!!!! Only keep data for links here !!!!
        MPI_Bcast(numlinks,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        if(my_rank != 0)
        {
            *ids = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
            *locs = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
            *numsteps = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
            data = (double***) malloc(*numlinks*sizeof(double**));
        }
        MPI_Bcast(*ids,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        MPI_Bcast(*locs,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        MPI_Bcast(*numsteps,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);

        for(i=0;i<*numlinks;i++)
        {
            if(my_rank != 0)
                data[i] = (double**) malloc((*numsteps)[i]*sizeof(double*));
            for(j=0;j<(*numsteps)[i];j++)
            {
                if(my_rank != 0)	data[i][j] = (double*) malloc(2*sizeof(double));
                MPI_Bcast(data[i][j],2,MPI_DOUBLE,0,MPI_COMM_WORLD);
            }
        }
    */

    MPI_Bcast(numlinks, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    if (my_rank != 0)
    {
        *numsteps = NULL;
        *ids = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
        *locs = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
    }
    MPI_Bcast(*ids, *numlinks, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(*locs, *numlinks, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    return data;
}



//Output functions ****************************************************************************
int Output_Linkid(unsigned int link_id, double t, VEC y_i, VEC global_params, VEC params, int state, void* user)
{
    CustomParams* forecastparams = (CustomParams*)user;
    return forecastparams->ID;
}

int Output_Timestamp(unsigned int link_id, double t, VEC y_i, VEC global_params, VEC params, int state, void* user)
{
    CustomParams* forecastparams = (CustomParams*)user;
    return (int)(round(t * 60.0 + forecastparams->offset - forecastparams->simulation_time_with_data) + 0.1);
}


//Custom parameters for forecasting ***********************************************************
void Init_Output_User_forecastparams(AsynchSolver* asynch)
{
    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
    Link* sys = asynch->sys;

    for (i = 0; i < my_N; i++)
        sys[my_sys[i]].output_user = malloc(sizeof(CustomParams));
}

void Free_Output_User_forecastparams(AsynchSolver* asynch)
{
    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
    Link* sys = asynch->sys;

    for (i = 0; i < my_N; i++)
    {
        free(sys[my_sys[i]].output_user);
        sys[my_sys[i]].output_user = NULL;
    }
}

void Set_Output_User_forecastparams(AsynchSolver* asynch, unsigned int offset, unsigned int time_with_data)
{
    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
    Link* sys = asynch->sys;
    CustomParams* forecastparams;

    for (i = 0; i < my_N; i++)
    {
        forecastparams = (CustomParams*)sys[my_sys[i]].output_user;
        forecastparams->ID = sys[my_sys[i]].ID;
        forecastparams->offset = offset;
        forecastparams->simulation_time_with_data = time_with_data;
    }
}

void Init_Output_PeakflowUser_Offset(AsynchSolver* asynch)
{
    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
    Link* sys = asynch->sys;

    for (i = 0; i < my_N; i++)
        sys[my_sys[i]].peakoutput_user = malloc(sizeof(unsigned int));
}

void Free_Output_PeakflowUser_Offset(AsynchSolver* asynch)
{
    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
    Link* sys = asynch->sys;

    for (i = 0; i < my_N; i++)
    {
        free(sys[my_sys[i]].peakoutput_user);
        sys[my_sys[i]].peakoutput_user = NULL;
    }
}

void Set_Output_PeakflowUser_Offset(AsynchSolver* asynch, unsigned int offset)
{
    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
    Link* sys = asynch->sys;

    for (i = 0; i < my_N; i++)
        *(unsigned int*)(sys[my_sys[i]].peakoutput_user) = offset;
}

