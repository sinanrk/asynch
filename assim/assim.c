#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

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

#include "asynch_interface.h"
#include "assim_ls_methods.h"
#include "assim_models.h"
#include "forecaster_methods.h"
#include "assim.h"
#include "optparse.h"


// Global variables
bool verbose = false;
int my_rank;
int np;

////Output functions
//int Output_Linkid(double t, VEC y_i, VEC global_params, VEC params, int state, void* user);
//void Set_Output_User_LinkID(asynchsolver* asynch);
//
//int Output_Timestamp(double t, VEC y_i, VEC global_params, VEC params, int state, void* user);
//void Set_Output_User_Timestamp(asynchsolver* asynch);


//Print to stdout only for process of rank 0
int print_out(const char* format, ...)
{
    int res = 0;
    if (my_rank == 0)
    {
        va_list args;
        va_start(args, format);
        res = vprintf(format, args);
        va_end(args);
    }

    return res;
}


//Print to stderr only for process of rank 0
int print_err(const char* format, ...)
{
    int res = 0;
    if (my_rank == 0)
    {
        va_list args;
        va_start(args, format);
        res = vfprintf(stderr, format, args);
        va_end(args);
    }

    return res;
}


//Make sure we finalize MPI
void asynch_onexit(void)
{
    int flag;
    MPI_Finalized(&flag);
    if (!flag)
        MPI_Finalize();
}


int main(int argc, char* argv[])
{
    int res;

    //Initialize MPI stuff
    res = MPI_Init(&argc, &argv);
    if (res == MPI_SUCCESS)
        atexit(asynch_onexit);
    else
    {
        print_err("Failed to initialize MPI");
        exit(EXIT_FAILURE);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //PETSC
    PetscInitialize(&argc, &argv, NULL, NULL);

    //Command line options
    bool debug = false;
    bool help = false;
    bool version = false;

    //Parse command line
    struct optparse options;
    optparse_init(&options, argv);
    struct optparse_long longopts[] = {
        { "debug", 'd', OPTPARSE_NONE },
        { "help", 'h', OPTPARSE_NONE },
        { "version", 'v', OPTPARSE_NONE },
        { "verbose", 'w', OPTPARSE_NONE },
        { 0 }
    };
    int option;
    while ((option = optparse_long(&options, longopts, NULL)) != -1) {
        switch (option) {
        case 'd':
            debug = true;
            break;
        case 'h':
            help = true;
            break;
        case 'v':
            version = true;
            break;
        case 'w':
            verbose = true;
            break;
        case '?':
            print_err("%s: %s\n", argv[0], options.errmsg);
            exit(EXIT_FAILURE);
        }
    }

    if (version) print_out("This is %s\n", PACKAGE_STRING);
    if (help)
    {
        print_out("Usage: asynch <global file>\n", PACKAGE_STRING);
        print_out(
            "  -d [--debug]   : Wait for the user input at the begining of the program (useful" \
            "                   for attaching a debugger)\n" \
            "  -w [--verbose] : Print debugging information to stdout\n" \
            "  -v [--version] : Print the current version of ASYNCH\n");
        exit(EXIT_SUCCESS);
    }
    if (version || help) exit(EXIT_SUCCESS);

    //Parse remaining arguments
    char *global_filename = optparse_arg(&options);
    if (global_filename == NULL && my_rank == 0)
    {
        print_err("Command line parameter required:  A universal variable file (.gbl).\n");
        exit(EXIT_FAILURE);
    }

    char *assim_filename = optparse_arg(&options);
    if (assim_filename == NULL && my_rank == 0)
    {
        print_err("Command line parameter required:  An assim file (.das).\n");
        exit(EXIT_FAILURE);
    }

    if (debug)
    {
        //Disable stdout buffering
        setvbuf(stdout, NULL, _IONBF, 0);

        //When the program first starts to execute, at the very beginning of our program, we 
        //ask the user to type some sort of input to simply stall the application until start your
        //"Attach to Process" and you can attach to all the different threads in your program.
        if (my_rank == 0)
        {
            printf("You may now attach the debugger then press enter.\n");
            //fflush(stdout);
            int ch = getchar();
        }

        MPI_Barrier(MPI_COMM_WORLD); // All threads will wait here until you give thread 0 an input
    }

    //Declare variables
    unsigned int i, j/*, k,l,m,n,*/;
    //RKMethod** AllMethods;
    //char additional[16];	//For output filename
    //srand(time(NULL));	//!!!! Is this needed? !!!!

    //Init asynch object and the river network
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
    Asynch_Parse_GBL(asynch, global_filename);
    if (my_rank == 0)	printf("Loading network...\n");
    Asynch_Load_Network(asynch);

    //Read data assimilation file
    //Create assim
    AssimData assim;
    InitAssimData(&assim, assim_filename, asynch);

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
    FindUpstreamLinks(asynch, &assim, problem_dim, true, assim.obs_time_step, assim.num_steps, assim.obs_locs, assim.num_obs);

    print_out("Partitioning network...\n");
    Asynch_Partition_Network(asynch);
    CleanUpstreamLinks(asynch);
    print_out("Loading parameters...\n");
    Asynch_Load_Network_Parameters(asynch, 0);
    print_out("Reading dam and reservoir data...\n");
    Asynch_Load_Dams(asynch);
    print_out("Setting up numerical error data...\n");
    Asynch_Load_Numerical_Error_Data(asynch);
    print_out("Initializing model...\n");
    Asynch_Initialize_Model(asynch);
    Setup_Errors(asynch, problem_dim);
    print_out("Loading initial conditions...\n");
    Asynch_Load_Initial_Conditions(asynch);
    print_out("Loading forcings...\n");
    Asynch_Load_Forcings(asynch);
    print_out("Loading output data information...\n");
    Asynch_Load_Save_Lists(asynch);
    print_out("Finalizing network...\n");
    Asynch_Finalize_Network(asynch);
    print_out("Calculating initial step sizes...\n");
    Asynch_Calculate_Step_Sizes(asynch);

    // No output needed, since we are not forecasting

    ////Setup output for link id, if needed
    //int setup_id = Asynch_Check_Output(asynch, "LinkID");
    //int setup_timestamp = Asynch_Check_Output(asynch, "Timestamp");
    //if (setup_id || setup_timestamp)
    //{
    //    if (my_rank == 0)	printf("Error: forecaster needs LinkID (%i), Timestamp (%i).\n", setup_id, setup_timestamp);
    //    MPI_Abort(MPI_COMM_WORLD, 1);
    //}

    ////Setup output for link id, if needed
    //int setup = Asynch_Check_Output(asynch, "LinkID");
    //if (setup != -1)
    //{
    //    Set_Output_User_LinkID(asynch);
    //    Asynch_Set_Output(asynch, "LinkID", ASYNCH_INT, &Output_Linkid, NULL, 0);
    //}

    ////Setup output for timestamp, if needed
    //setup = Asynch_Check_Output(asynch, "Timestamp");
    //if (setup != -1)
    //{
    //    Asynch_Set_Output(asynch, "Timestamp", ASYNCH_INT, &Output_Timestamp, NULL, 0);
    //}

    ////Prepare output files
    //Asynch_Prepare_Temp_Files(asynch);
    //Asynch_Write_Current_Step(asynch);
    //Asynch_Prepare_Peakflow_Output(asynch);
    //Asynch_Prepare_Output(asynch);

    //Pull data from asynch
    unsigned int my_N = asynch->my_N, N = asynch->N, *my_sys = asynch->my_sys, **id_to_loc = asynch->id_to_loc, num_forcings = asynch->globals->num_forcings;
    int *assignments = asynch->assignments;
    Link* sys = asynch->sys;
    short int *getting = asynch->getting;
    GlobalVars *globals = asynch->globals;
    Model* custom_model = asynch->custom_model;

    //Set print_time to t_0
    //Asynch_Reset_Temp_Files(asynch,sys[my_sys[0]]->last_t);

    ////Reserve space for backups
    //VEC* backup = (VEC*)malloc(N * sizeof(VEC));	//!!!! Same as background x_b? !!!!
    //for (i = 0; i < N; i++)
    //{
    //    if (assignments[i] == my_rank || getting[i])
    //        backup[i] = v_get(sys[i].dim);
    //    else
    //        backup[i] = v_get(0);
    //}

    //Initialize choices
    unsigned int num_obs = assim.num_obs, *obs_locs = assim.obs_locs, num_steps = assim.num_steps;
    double t_b = 0.0, obs_time_step = assim.obs_time_step;
    double t_f = Asynch_Get_Total_Simulation_Duration(asynch);
    unsigned int allstates = assim_dim * N;
    //double x_b[allstates];

    // Allocate background
    double  *x_b = calloc(allstates, sizeof(double));
    for (i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank)
        {
            for (j = 0; j < assim_dim; j++)
                x_b[i*assim_dim + j] = sys[i].list->tail->y_approx.ve[j];	//!!!! Need to be able to specify which states are used !!!!
        }
        //MPI_Bcast(&(x_b[i*assim_dim]), assim_dim, MPI_DOUBLE, assignments[i], MPI_COMM_WORLD);
    }

    int mpi_res = MPI_Allreduce(MPI_IN_PLACE, x_b, allstates, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //Other initializations
    //unsigned int numstepstotake;
    //unsigned int iterations = (unsigned int)round((t_f - t_b) / inc);

    // Vector of observations (single step)
    double* d = calloc(assim.num_obs, sizeof(double));

    // Vector of observations (multiple steps)
    double* d_full = calloc(num_obs * num_steps, sizeof(double));

    //Values used to start asynch solver in tao solvers
    double* x_start = calloc(allstates, sizeof(double));	//Values used to start asynch solver in tao solvers

    //Call model specific data assimilation routines
    //For Model 254
    //unsigned int allstates_needed = Setup_Fitting_Data_Model254(asynch,obs_locs,num_obs);
    //For Model 254 trim, q
    Setup_Fitting_Data_Model254_q(asynch, obs_locs, num_obs);
    //For Model 254 trim, q and s_p
    //Setup_Fitting_Data_Model254_qsp(asynch,obs_locs,num_obs);
    //For Model 254 trim, q and s_t
    //Setup_Fitting_Data_Model254_qst(asynch,obs_locs,num_obs);

    //Find locations unaffected by gauges
    unsigned int *vareq_shift, *inv_vareq_shift;
    unsigned int allstates_needed = BuildStateShift(asynch, allstates, obs_locs, num_obs, &vareq_shift, &inv_vareq_shift);


    printf("allstates_needed: %u allstates: %u\n", allstates_needed, allstates);

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

    //Prep PetSC
    if (my_rank == 0)
        printf("\nPrepping PetSc...\n");
    Workspace ws;

    //For linear least squares
    int *HM_col_indices = NULL;
    if (my_rank == 0)
    {
        VecCreateSeq(MPI_COMM_SELF, allstates_needed, &ws.rhs);
        VecCreateSeq(MPI_COMM_SELF, allstates_needed, &ws.x);
        MatCreateSeqDense(MPI_COMM_SELF, num_obs*num_steps, allstates_needed, NULL, &ws.HM);
        MatCreateSeqDense(MPI_COMM_SELF, allstates_needed, allstates_needed, NULL, &ws.HTH);
        MatCreateSeqDense(MPI_COMM_SELF, allstates_needed, num_obs*num_steps, NULL, &ws.HMTR);
        HM_col_indices = (int*)malloc(allstates_needed * sizeof(int));
        for (i = 0; i < allstates_needed; i++)
            HM_col_indices[i] = i;
        KSPCreate(MPI_COMM_SELF, &ws.ksp);
        KSPSetOperators(ws.ksp, ws.HTH, ws.HTH);
        //KSPSetTolerances(ws.ksp, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
        KSPSetFromOptions(ws.ksp);	// This is used to override the solver setting from the command line
    }

    int* d_indices = (int*)malloc(num_steps*num_obs * sizeof(int));
    for (i = 0; i < num_steps*num_obs; i++)
        d_indices[i] = i;

    //Transfer upstreams areas to all procs
    double* inv_upareas = (double*)calloc(N, sizeof(double));
    UpstreamData* updata;
    for (i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank)
            inv_upareas[i] = 1.0 / (sys[i].params.ve[globals->area_idx] * 1e3);
    }

    MPI_Allreduce(MPI_IN_PLACE, inv_upareas, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //Links needed for fitting
    bool *links_needed = (bool*)calloc(N, sizeof(bool));
    for (i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == my_rank)
        {
            Link *current = &sys[obs_locs[i]];
            updata = (UpstreamData*)(current->user);
            links_needed[current->location] = true;
            //for(j=0;j<current->num_parents;j++)
            //{
            //	for(k=0;k<updata->num_upstreams[j];k++)
            //		my_links_needed[updata->upstreams[j][k]] = 1;
            //}

            for (j = 0; j < updata->num_upstreams; j++)
            {
                links_needed[updata->upstreams[j]->location] = true;
            }
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, links_needed, N, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

    //Build weight matrices
    //!!!! Assuming only q is changing !!!!
    unsigned int curr_idx = 0;
    if (my_rank == 0)
    {
        VecCreateSeq(MPI_COMM_SELF, allstates_needed, &ws.B);
        for (i = 0; i < N; i++)
        {
            if (links_needed[i])
                //VecSetValue(ws.B, curr_idx++, inv_upareas[i], INSERT_VALUES);
                VecSetValue(ws.B, curr_idx++, 1.0, INSERT_VALUES);
        }
        VecAssemblyBegin(ws.B);
        VecAssemblyEnd(ws.B);

        VecCreateSeq(MPI_COMM_SELF, num_obs*num_steps, &ws.R);
        for (i = 0; i < num_obs; i++)
        {
            for (j = 0; j < num_steps; j++)
                //VecSetValue(ws.R, j*num_obs + i, inv_upareas[obs_locs[i]] * 10.0, INSERT_VALUES);
                VecSetValue(ws.R, j*num_obs + i, 1.0, INSERT_VALUES);
        }
        VecAssemblyBegin(ws.R);
        VecAssemblyEnd(ws.R);

        if (verbose)
        {
            printf("Weighting Matrix B (diagonal)\n");
            VecView(ws.B, PETSC_VIEWER_STDOUT_SELF);
            printf("Weighting Matrix R (diagonal)\n");
            VecView(ws.R, PETSC_VIEWER_STDOUT_SELF);
        }
    }
    free(links_needed);

    ws.HM_buffer = (double*)calloc(allstates_needed, sizeof(double));
    ws.HM_col_indices = HM_col_indices;
    ws.d_indices = d_indices;
    ws.d_full = d_full;
    ws.x_start = x_start;
    ws.asynch = asynch;
    ws.problem_dim = problem_dim;
    ws.assim_dim = assim_dim;
    ws.allstates = allstates;
    ws.allstates_needed = allstates_needed;
    ws.vareq_shift = vareq_shift;
    ws.inv_vareq_shift = inv_vareq_shift;
    ws.obs_time_step = obs_time_step;
    ws.num_steps = assim.num_steps;
    ws.obs_locs = obs_locs;
    ws.num_obs = num_obs;
    ws.t_b = t_b;
    ws.x_b = x_b;

    //Print out some information
    unsigned int my_eqs = 0, total_eqs;
    for (i = 0; i < my_N; i++)	my_eqs += sys[my_sys[i]].dim;
    MPI_Reduce(&my_eqs, &total_eqs, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

    printf("[%i]: Good to go with %u links (%u eqs).\n", my_rank, my_N, my_eqs);
    if (my_rank == 0)
    {
        ASYNCH_SLEEP(1);
        printf("\nNetwork has a total of %u links and %u equations.\n\n", N, total_eqs);
        printf("Making calculations...\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //// Backup the solution
    //k = 0;
    //for (i = 0; i < N; i++)
    //    if (backup[i].dim > 0)
    //        v_copy(asynch->sys[i].list->tail->y_approx, backup[i]);

    //double simulation_time_with_data = asynch->forcings[forecast_idx]->file_time * forecaster->num_rainsteps;
    double simulation_time_with_data = (assim.num_steps - 1) * assim.obs_time_step;
    unsigned int simulation_time_with_data_secs = (int)(simulation_time_with_data + 1e-3) * 60;

    //Setup temp files
    //Set_Output_User_forecastparams(asynch, forecast_time_unix, simulation_time_with_data_secs);	//!!!! Should this be done here at all? !!!!
    //Set_Output_PeakflowUser_Offset(asynch,forecast_time_unix);
    //Asynch_Set_Total_Simulation_Time(asynch, forecast_window + simulation_time_with_data);
    //Asynch_Prepare_Temp_Files(asynch);

    //if (my_rank == 0)
    //{
    //    time_t now = time(NULL);
    //    struct tm* now_info = localtime(&now);
    //    printf("\n\nPass %u\n", k);
    //    printf("Current time is %s", asctime(now_info));
    //}

    //Clear buffers
    Flush_TransData(asynch->my_data);

    //Make some initializations
    //first_file = last_file;
    //last_file = last_file + (unsigned int) asynch->forcings[forecast_idx]->file_time * 60 * num_rainsteps;
    //forecast_time_unix += change_time;
    //unsigned int forecast_idx = 0;
    //unsigned int forecast_time_unix = (unsigned int) asynch->forcings[forecast_idx]->first_file;
    //unsigned int background_time_unix = forecast_time_unix - (assim.num_steps - 1) * (unsigned int)(assim.obs_time_step * 60.0 + 1.0e-3);
    //nextforcingtime = forecast_time_unix - 60 * (unsigned int)asynch->forcings[forecast_idx]->file_time;	//This is the actual timestamp of the last needed forcing data. This will be downloaded (unlike last_file)


    unsigned int forcing_idx_rain = 0;
    unsigned int forcing_idx_tep = 1;

    unsigned int begin_assim_window = asynch->forcings[forcing_idx_rain].first_file;
    unsigned int end_assim_window = asynch->forcings[forcing_idx_rain].first_file + (unsigned int)(assim.num_steps * assim.obs_time_step * 60.0);


    //Reset each link
    //Asynch_Set_System_State(asynch, 0.0, backup);
    //Set_Output_User_forecastparams(asynch, forecast_time_unix, simulation_time_with_data_secs);
    //Set_Output_PeakflowUser_Offset(asynch,forecast_time_unix);
    //Asynch_Write_Current_Step(asynch);
    Asynch_Set_Forcing_State(asynch, forcing_idx_rain, 0.0, begin_assim_window, end_assim_window);

    for (i = 0; i < asynch->globals->num_forcings; i++)	//Set any other database forcings to begin at first_file
    {
        if (asynch->forcings[i].flag == 3)
            Asynch_Set_Forcing_State(asynch, i, 0.0, begin_assim_window, end_assim_window);
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

    double start = MPI_Wtime();





    //Start the analysis




    //time_t q_start, q_stop;
    //int *assignments = asynch->assignments;
    //short int *getting = asynch->getting;
    //unsigned int **id_to_loc = asynch->id_to_loc, N = asynch->N, num_obs = assim.num_obs;
    //Link **sys = asynch->sys;
    //UnivVars* globals = asynch->globals;
    //unsigned int i, j, l;
    //unsigned int steps_to_use = assim.steps_to_use;
    //unsigned int assim_window_unix = assim.num_steps * (int)(assim.obs_time_step + 1e-3);
    //double t_b = 0.0, inc = assim.inc;
    //double forecast_window = forecaster->forecast_window;
    //double assim_window = assim.num_steps * assim.obs_time_step;
    //double *d_full = ws->d_full, *x_start = ws->x_start, *x_b = ws->x_b;
    //unsigned int allstates = ws->allstates;
    unsigned int least_squares_iters = assim.least_squares_iters;
    double *analysis = (double*)malloc(allstates * sizeof(double));	//!!!! Should be removed !!!!
    //model* custom_model = asynch->custom_model;
    //double q[steps_to_use*num_obs];
    double *q = (double*)malloc(num_steps*num_obs * sizeof(double));


    //Set the forecast window
    Asynch_Set_Total_Simulation_Duration(asynch, num_steps * obs_time_step);	//!!!! Is this needed? !!!!

    //Get the new observations
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (my_rank == 0)
            printf("Downloading observations...\n");

        double start = MPI_Wtime();

        while (GetObservationsData(&assim, id_to_loc, N, begin_assim_window, d_full) == -1)
        {
            if (my_rank == 0)	printf("Error downloading observations. Retrying...\n");
            ASYNCH_SLEEP(5);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        double stop = MPI_Wtime();

        if (my_rank == 0)
            printf("Time to get new discharges: %.0f\n", stop - start);
    }

    if (verbose && my_rank == 0)
    {
        printf("d_full\n");
        Print_VECTOR(d_full, num_steps*num_obs);
        printf("\n");
    }

    ////Scale by upstream area
    //printf("Scaling init discharges by upstream area...\n");
    //AdjustDischarges(asynch, obs_locs, d_full, num_obs, assim_dim, x_b);

    //Set the initial guess from background
    memcpy(x_start, x_b, allstates * sizeof(double));
    //Copy in states that won't be affected by the optimization solver
    memcpy(analysis, x_b, allstates * sizeof(double));

    if (verbose && my_rank == 0)
    {
        //printf("Going in, minimizer...\n");
        //Print_VECTOR(minimizer,allstates_needed);
        //printf("Going in, x_start...\n");
        //Print_VECTOR(x_start, allstates);

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
    bool try_again = false;
    do
    {
        int iterations = 0;
        double error, prev_error = -1.0;
        for (j = 0; j < least_squares_iters; j++)
            //while(diff > 1e-2)
        {
            iterations++;
            LinearLeastSquares(&ws, q);
            error = compute_diff(d_full, q, num_steps*num_obs);
            if (prev_error >= 0.0)
            {
                double diff = prev_error - error;
                if (error > prev_error)
                {
                    if (my_rank == 0)
                    {
                        printf("!!!! LS error got worse. Breaking... !!!!\n");
                        printf("Errors are %f and %f\n", error, prev_error);
                    }

                    //Go back to previous solution
                    for (i = 0; i < allstates; i++)
                        x_start[i] = analysis[i];

                    break;
                }
                if (my_rank == 0)	printf("Difference is %f (%f vs %f)\n", diff, error, prev_error);
            }

            prev_error = error;
            for (i = 0; i < allstates; i++)
                analysis[i] = x_start[i];
        }
        if (my_rank == 0)
            printf("Total iterations = %i\n", iterations);

        //if(error > 30.0 && !try_again)	//Check if the numerical scheme is having convergence issues
        //if(error > 10.0 && !try_again)
        if (!try_again)
        {
            //try_again = 1;
            try_again = ReduceBadDischargeValues(sys, assignments, N, d_full, q, num_steps, ws.obs_locs, num_obs, x_start, assim_dim, 1.0);	//!!!! Not sure what to use for the limit... !!!!
        }
        else
            try_again = false;
    } while (try_again);


    //Copy x_start to analysis  !!!! This shouldn't happen. Try a pointer dance. Or maybe something better... !!!!
    //for(i=0;i<allstates;i++)	analysis[i] = x_start[i];

    if (verbose && my_rank == 0)
    {
        //if(k == 10)
        {
            printf("x_b\n");
            Print_VECTOR(x_b, allstates);
            printf("\n");

            printf("analysis [%d - %d]\n", begin_assim_window, end_assim_window);
            Print_VECTOR(analysis, allstates);
        }
    }

    free(analysis);
    free(q);
















    double stop = MPI_Wtime();
    print_out("\nTime for calculations: %f. All done!\n", stop - start);

    //Prepare snapshots
    for (i = 0; i < N; i++)
    {
        Link *current = &asynch->sys[i];
        if (current->list != NULL)
            for (j = 0; j < problem_dim; j++)
                current->list->head->y_approx.ve[j] = x_start[i * problem_dim + j];
    }

    //Make a snaphsot
    print_out("Making snapshot\n");
    Asynch_Take_System_Snapshot(asynch, NULL);

    //Clean up
    print_out("Cleaning up\n");
    free(d);
    free(d_full);
    free(x_start);
    free(vareq_shift);
    free(inv_vareq_shift);

    if (my_rank == 0)
    {
        MatDestroy(&ws.HM);
        MatDestroy(&ws.HTH);
        VecDestroy(&ws.rhs);
        VecDestroy(&ws.x);
        VecDestroy(&ws.R);
        VecDestroy(&ws.B);
        MatDestroy(&ws.HMTR);
        KSPDestroy(&ws.ksp);
    }
    free(HM_col_indices);
    free(d_indices);

    free(inv_upareas);
    FreeAssimData(&assim);
    //Free_ForecastData(&forecaster);
    //Free_Output_PeakflowUser_Offset(asynch);
    //Free_Output_User_forecastparams(asynch);

    //Petsc clean up
    PetscFinalize();

    //Asynch clean up
    FreeUpstreamLinks(asynch);
    Asynch_Delete_Temporary_Files(asynch);
    Asynch_Free(asynch);

    return res;
}


//void Make_Assimilated_Forecasts(asynchsolver* asynch, unsigned int background_time_unix, double simulation_time_with_data, VEC* backup, AppCtxFull* user, ForecastData* forecaster, AssimData* assim, unsigned int assim_dim, unsigned int forecast_idx)
//{
//    time_t q_start, q_stop;
//    int *assignments = asynch->assignments;
//    short int *getting = asynch->getting;
//    unsigned int **id_to_loc = asynch->id_to_loc, N = asynch->N, num_obs = assim->num_obs;
//    Link **sys = asynch->sys;
//    UnivVars* globals = asynch->globals;
//    unsigned int i, j, l, steps_to_use = assim->steps_to_use, assim_window_unix = assim->steps_to_use * (int)(assim->inc + 1e-3);
//    double t_b = 0.0, inc = assim->inc, forecast_window = forecaster->forecast_window, assim_window = assim->steps_to_use * assim->inc;
//    double *d_full = user->d_full, *x_start = user->x_start, *x_b = user->x_b;
//    unsigned int allstates = user->allstates, least_squares_iters = assim->least_squares_iters;
//    double *analysis = (double*)malloc(allstates * sizeof(double));	//!!!! Should be removed !!!!
//    model* custom_model = asynch->custom_model;
//    //double q[steps_to_use*num_obs];
//    double *q = (double*)malloc(steps_to_use*num_obs * sizeof(double));
//
//    //for(k=0;k<iterations;k++)
//    {
//        /*
//                if(my_rank == 0)
//                {
//                    printf("\n\n*************************************\n");
//                    printf("Iteration %u/%u. Background = %e.\n",k,iterations,t_b);
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//
//                steps_to_use = (k < max_steps_to_use) ? k+1 : max_steps_to_use;
//        */
//        //Set the forecast window
//        Asynch_Set_Total_Simulation_Time(asynch, forecast_window + steps_to_use*inc);	//!!!! Is this needed? !!!!
//
//        //Get the new observations
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (my_rank == 0)
//            printf("Downloading observations...\n");
//        time(&q_start);
//
//        while (GetObservationsData(assim, id_to_loc, N, background_time_unix, d_full) == -1)
//        {
//            if (my_rank == 0)	printf("Error downloading observations. Retrying...\n");
//            ASYNCH_SLEEP(5);
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//
//        /*
//                !!!! Need observations !!!!
//                GetObseravations(assim,background_time_unix + assim_window_unix,d);	!!!! Only gets newest. Need more for first time through !!!!
//                //FindAllDischarges(truesolution,0.0 + (k)*inc,num_obs,numsteps,d);
//                for(j=0;j<(steps_to_use-1)*num_obs;j++)
//                    d_full[j] = d_full[j+num_obs];
//                for(j=0;j<num_obs;j++)
//                    d_full[(steps_to_use-1)*num_obs + j] = d[j];
//        */
//        MPI_Barrier(MPI_COMM_WORLD);
//        time(&q_stop);
//        if (my_rank == 0)
//            printf("Time to get new discharges: %.0f\n", difftime(q_stop, q_start));
//
//
//        if (verbose && my_rank == 0)
//        {
//            printf("d_full\n");
//            Print_VECTOR(d_full, steps_to_use*num_obs);
//            printf("\n");
//        }
//
//
//        /*
//        //if(k == 0)	//!!!! Move this outside the loop. Can I just change minimizer? !!!!
//        {
//            AdjustDischarges(asynch,assim->obs_locs,d_full,num_obs,x_b,allstates,assim_dim);
//        }
//        */
//
//
//        //Set the initial guess
//        for (i = 0; i < N; i++)	//Copy in states that won't be affected by the optimization solver
//        {
//            for (j = 0; j < assim_dim; j++)
//                x_start[i*assim_dim + j] = analysis[i*assim_dim + j] = x_b[i*assim_dim + j];	//!!!! Do I need all these? !!!!
//        }
//
//
//        if (verbose && my_rank == 0)
//        {
//            //printf("Going in, minimizer...\n");
//            //Print_VECTOR(minimizer,allstates_needed);
//            printf("Going in, x_start...\n");
//            Print_VECTOR(x_start, allstates);
//
//            /*
//            double* tempy;
//            printf("Guess:\n");
//            VecGetArray(Guess,&tempy);
//            for(i=0;i<allstates_needed;i++)
//            printf("%.15e ",tempy[i]);
//            printf("\n");
//            VecRestoreArray(Guess,&tempy);
//            */
//        }
//
//        //Calculate the analysis
//        int try_again = 0;
//        do
//        {
//            int iterations = 0;
//            double diff = 10.0, error, prev_error = -1.0;
//            for (j = 0; j < least_squares_iters; j++)
//                //while(diff > 1e-2)
//            {
//                iterations++;
//                LinearLeastSquares(user, q);
//                error = compute_diff(d_full, q, steps_to_use*num_obs);
//                if (prev_error >= 0.0)
//                {
//                    diff = prev_error - error;
//                    if (error > prev_error)
//                    {
//                        if (my_rank == 0)
//                        {
//                            printf("!!!! LS error got worse. Breaking... !!!!\n");
//                            printf("Errors are %f and %f\n", error, prev_error);
//                        }
//
//                        //Go back to previous solution
//                        for (i = 0; i < allstates; i++)	x_start[i] = analysis[i];
//
//                        break;
//                    }
//                }
//                if (my_rank == 0)	printf("Difference is %f (%f vs %f)\n", diff, error, prev_error);
//                prev_error = error;
//                for (i = 0; i < allstates; i++)	analysis[i] = x_start[i];
//            }
//            if (my_rank == 0)
//                printf("Total iterations = %i\n", iterations);
//
//            //if(error > 30.0 && !try_again)	//Check if the numerical scheme is having convergence issues
//            //if(error > 10.0 && !try_again)
//            if (!try_again)
//            {
//                //try_again = 1;
//                try_again = ReduceBadDischargeValues(sys, assignments, N, d_full, q, steps_to_use, user->obs_locs, num_obs, x_start, assim_dim, 1.0);	//!!!! Not sure what to use for the limit... !!!!
//            }
//            else	try_again = 0;
//        } while (try_again);
//
//
//        //Copy x_start to analysis  !!!! This shouldn't happen. Try a pointer dance. Or maybe something better... !!!!
//        //for(i=0;i<allstates;i++)	analysis[i] = x_start[i];
//
//        if (verbose && my_rank == 0)
//        {
//            //if(k == 10)
//            {
//                printf("x_b\n");
//                Print_VECTOR(x_b, allstates);
//                printf("\n");
//
//                printf("analysis t_b = %e last assim time = %e\n", t_b, t_b + 5.0*(steps_to_use - 1));
//                Print_VECTOR(analysis, allstates);
//                //getchar();
//            }
//        }
//
//        //Switch to model without variational eqs
//        //For model 254
//        //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Model_254,&Precalculations_Assim_254,&ReadInitData_Assim_254);
//        //For model 254 trim
//        Asynch_Custom_Model(asynch, &SetParamSizes_Assim_254, &ConvertParams_Assim_254, &InitRoutines_Model_252, &Precalculations_Assim_254, &ReadInitData_Assim_254_q);	//!!!! I think this is ok... !!!!
//        Asynch_Initialize_Model(asynch);
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        time(&q_start);
//
//        //Advance to get new background
//        ResetSysLS(sys, N, globals, t_b, analysis, assim_dim, asynch->globals->num_forcings, asynch->my_data);
//        for (i = 0; i < N; i++)				//!!!! Put this into ResetSysLS? !!!!
//        {
//            if (assignments[i] == my_rank || getting[i])
//            {
//                Link *current = sys[i];
//                custom_model->InitializeEqs(globals->global_params, current->params, NULL, 0, current->list->head->y_approx, globals->type, current->diff_start, current->no_ini_start, current->user, NULL);
//            }
//        }
//
//        {
//            Asynch_Set_Total_Simulation_Time(asynch, simulation_time_with_data);
//            //globals->maxtime = t_b + inc;
//            Asynch_Advance(asynch, 1);
//
//            //Extract background !!!! Isn't there a routine for this? (there needs to be...) !!!!
//            for (j = 0; j < N; j++)
//            {
//                if (assignments[j] == my_rank)
//                {
//                    for (l = 0; l < assim_dim; l++)
//                        x_b[j*assim_dim + l] = backup[j].ve[l] = sys[j]->list->tail->y_approx.ve[l];
//                }
//                MPI_Bcast(&(x_b[j*assim_dim]), assim_dim, MPI_DOUBLE, assignments[j], MPI_COMM_WORLD);
//            }
//        }
//
//
//        //Make forecasts
//        Asynch_Deactivate_Forcing(asynch, forecast_idx);
//        Asynch_Set_Total_Simulation_Time(asynch, forecast_window + simulation_time_with_data);
//        Asynch_Advance(asynch, 1);
//        Asynch_Activate_Forcing(asynch, forecast_idx);
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        time(&q_stop);
//        if (my_rank == 0)	printf("Time for forecast: %.0f\n", difftime(q_stop, q_start));
//
//        //Switch back to model with variational eqs
//        //For model 254
//        //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254,&Precalculations_Assim_254,&ReadInitData_Assim_254);
//        //For model 254 trim, q
//        Asynch_Custom_Model(asynch, &SetParamSizes_Assim_254, &ConvertParams_Assim_254, &InitRoutines_Assim_254_q, &Precalculations_Assim_254, &ReadInitData_Assim_254_q);
//        //Model 254, q and s_p
//        //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_qsp,&Precalculations_Assim_254,&ReadInitData_Assim_254_qsp);
//        //Model 254, q and s_t
//        //Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_qst,&Precalculations_Assim_254,&ReadInitData_Assim_254_qst);
//        Asynch_Initialize_Model(asynch);
//        /*
//                //Go to next time
//                t_b += inc;
//                user.t_b = t_b;
//        */
//
//        /*
//                //Reset the temporary files
//                if(my_rank == 0)	printf("Creating output file...\n");
//        MPI_Barrier(MPI_COMM_WORLD);
//        time(&q_start);
//                sprintf(additional,"%u",k);	//Add the iteration number to the output files
//                Asynch_Create_Output(asynch,additional);
//                Asynch_Reset_Temp_Files(asynch,t_b);
//        MPI_Barrier(MPI_COMM_WORLD);
//        time(&q_stop);
//        if(my_rank == 0)	printf("Time to create output: %.0f\n",difftime(q_stop,q_start));
//        */
//    }
//
//    free(analysis);
//    free(q);
//}

//This computes the least squares fit assuming the background and analysis difference is linear in the innovations. 
//HM is (num_obs*steps_to_use) X allstates_needed
//HM_els is 1 X allstates (i.e. a 1D array)
//HM_buffer is 1 X allstates_needed
int LinearLeastSquares(Workspace* ws, double* q)
{
    unsigned int i, j,/*k,m,*/n/*,l,counter*/;

    //Unpack ptr
    AsynchSolver* asynch = ws->asynch;
    unsigned int N = asynch->N;
    Link* sys = asynch->sys;
    GlobalVars* globals = asynch->globals;
    int* assignments = asynch->assignments;
    short int* getting = asynch->getting;
    //Mat *HM = ws->HM, *HTH = ws->HTH, *HMTR = ws->HMTR;
    //Vec *RHS = ws->RHS, d, *x = ws->x, *B = ws->B, *R = ws->R;
    //KSP *ksp = ws->ksp;
    Vec d;
    unsigned int *obs_locs = ws->obs_locs, assim_dim = ws->assim_dim;
    double obs_time_step = ws->obs_time_step;
    unsigned int problem_dim = ws->problem_dim, allstates = ws->allstates;
    //unsigned int steps_to_use = ws->steps_to_use, num_obs = ws->num_obs;
    int *HM_col_indices = ws->HM_col_indices, *d_indices = ws->d_indices;
    double t_b = ws->t_b;
    unsigned int allstates_needed = ws->allstates_needed;
    double /**RHS_els,*/*x_start = ws->x_start, *HM_buffer = ws->HM_buffer;
    //unsigned int max_or_steps = steps_to_use;
    Model* custom_model = asynch->custom_model;
    unsigned int *vareq_shift = ws->vareq_shift, *inv_vareq_shift = ws->inv_vareq_shift;

    double start = MPI_Wtime();

    unsigned int num_total_obs = ws->num_steps * ws->num_obs;

    //Build a vector structure for d !!!! Obviously, this needs to not happen... !!!!
    //if(max_or_steps > allstates_needed)	printf("[%i]: Error: max_or_steps > allstates_needed (%u > %u)\n",my_rank,max_or_steps,allstates_needed);
    VecCreateSeq(MPI_COMM_SELF, num_total_obs, &d);	//!!!! This needs to be fixed !!!!
    VecSet(d, 0.0);
    VecSetValues(d, num_total_obs, d_indices, ws->d_full, INSERT_VALUES);
    VecAssemblyBegin(d);
    VecAssemblyEnd(d);

    /*
    //!!!! Try to scale the init conditions !!!!
    AdjustDischarges(asynch,obs_locs,d_els,num_obs,x_start,allstates);
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
    for (i = 0; i < ws->num_steps; i++)	//!!!! Start at i=1? For i = 0, I don't think we need to set anything... !!!!		//HM here holds the values of M that are needed
    {
        globals->maxtime = t_b + i * obs_time_step;
        if (i)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            double start = MPI_Wtime();

            Asynch_Advance(asynch, 0);

            MPI_Barrier(MPI_COMM_WORLD);
            double stop = MPI_Wtime();

            if (my_rank == 0)
                printf("Time for advance to time %f: %.0f\n", globals->maxtime, stop - start);
        }


        //printf("ID = %u, t = %e\n",sys[obs_locs[0]]->ID,sys[obs_locs[0]]->last_t);
        //Print_Vector(sys[obs_locs[0]]->list->tail->y_approx);
        //printf("**********\n");

        //Build HM
        for (j = 0; j < ws->num_obs; j++)
        {
            Link *current = &sys[obs_locs[j]];	//!!!! Assumes only discharges !!!!
            int owner = assignments[obs_locs[j]];
            bool is_my_link = (owner == my_rank);

            UpstreamData *updata = (UpstreamData*)(current->user);
            memset(HM_buffer, 0, allstates_needed * sizeof(double));

            //From my link
            if (is_my_link)
            {
                //Pull out needed data
                for (n = 0; n < updata->num_fit_states; n++)
                {
                    printf("ID = %u | Loading %e (from %u) into spot %u\n",
                        current->ID,
                        current->list->tail->y_approx.ve[updata->fit_states[n]],
                        updata->fit_states[n],
                        vareq_shift[updata->fit_to_universal[n]]);
                    //ASYNCH_SLEEP(1);
                    HM_buffer[vareq_shift[updata->fit_to_universal[n]]] = current->list->tail->y_approx.ve[updata->fit_states[n]];
                }

                //Extract calculationed q's (Just needed for testing. Maybe...)
                q[i * ws->num_obs + j] = current->list->tail->y_approx.ve[0];
            }

            //MPI_Bcast(HM_buffer, allstates_needed, MPI_DOUBLE, owner, MPI_COMM_WORLD);	//!!!! Only proc 0 needs this !!!!
            if (my_rank == 0)
            {
                MPI_Reduce(MPI_IN_PLACE, HM_buffer, allstates_needed, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

                unsigned int row_idx = i * ws->num_obs + j;
                MatSetValues(ws->HM, 1, &row_idx, allstates_needed, HM_col_indices, HM_buffer, INSERT_VALUES);
            }
            else
                MPI_Reduce(HM_buffer, NULL, allstates_needed, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



            MPI_Bcast(&(q[i * ws->num_obs + j]), 1, MPI_DOUBLE, owner, MPI_COMM_WORLD);

            //printf("Got %u, %u\n",allstates_needed,updata->num_fit_states);
            //for(n=0;n<allstates_needed;n++)
            //	printf("%e ",HM_buffer[n]);
            //printf("\n");
            //char r = getchar();                
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
                MatSetValues(ws->HM,1,&i,allstates,cols_allstates_needed,HM_els,INSERT_VALUES);
        }
    */

    if (my_rank == 0)
    {
        //Assemble the HM matrix
        MatAssemblyBegin(ws->HM, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(ws->HM, MAT_FINAL_ASSEMBLY);

        if (verbose)
        {
            printf("Matrix HM\n");
            MatView(ws->HM, PETSC_VIEWER_STDOUT_SELF);
        }

        start = MPI_Wtime();

        //Calculate innovations
        double *buffer = NULL;
        VecGetArray(d, &buffer);
        for (i = 0; i < num_total_obs; i++)
            buffer[i] = buffer[i] - q[i];
        VecRestoreArray(d, &buffer);

        //Build the linear system \f$ A x = rhs \f$
        //HMTR is allstates_needed x (num_obs*max_or_steps)
        //HM is (num_obs*max_or_steps) x allstates_needed

        /// \f$ HMTR = H(y_0)^T R \f$
        /// HMTR is a temporary variable used for rhs and A computation
        MatTranspose(ws->HM, MAT_REUSE_MATRIX, &ws->HMTR);
        MatDiagonalScale(ws->HMTR, NULL, ws->R);

        /// \f$ A = B + H(y_0)^T R H(y_0) \f$
        MatMatMult(ws->HMTR, ws->HM, MAT_REUSE_MATRIX, PETSC_DEFAULT, &ws->HTH);
        MatDiagonalSet(ws->HTH, ws->B, ADD_VALUES);
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

        /// \f$ rhs = H(y_0)^T R \alpha(y_0^b) \f$
        MatMult(ws->HMTR, d, ws->rhs);
        MatAssemblyBegin(ws->HTH, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(ws->HTH, MAT_FINAL_ASSEMBLY);

        if (verbose)
        {
            printf("Matrix HTH\n");
            MatView(ws->HTH, PETSC_VIEWER_STDOUT_SELF);
        }

        //MPI_Barrier(MPI_COMM_WORLD);
        double stop = MPI_Wtime();
        if (my_rank == 0)
            printf("Time for matrix computations: %.0f\n", stop - start);

        //Compute analysis
    //MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        /// \f$ x = y_0 - y_0^b \f$
        KSPSetOperators(ws->ksp, ws->HTH, ws->HTH);     //Maybe not actually necessary
        KSPSolve(ws->ksp, ws->rhs, ws->x);
        KSPConvergedReason reason;
        KSPGetConvergedReason(ws->ksp, &reason);

        if (my_rank == 0)
            printf("Converged reason: %s\n", KSPConvergedReasons[reason]);

        //MPI_Barrier(MPI_COMM_WORLD);
        stop = MPI_Wtime();
        if (my_rank == 0)
            printf("Time for inversion: %.0f\n", stop - start);

        if (verbose)
        {
            printf("Solution x\n");
            VecView(ws->x, PETSC_VIEWER_STDOUT_SELF);
        }

        //Copy new solution to x_start
        VecGetArray(ws->x, &buffer);
        //for(i=0;i<num_above;i++)
        for (i = 0; i < allstates_needed; i++)	//!!!! I think this is right... !!!!
        {
            //printf("i = %u %u\n",i,inv_vareq_shift[i]);
            //ASYNCH_SLEEP(1);
            x_start[inv_vareq_shift[i]] += buffer[i];
            //x_start[above_gauges[i]*assim_dim] += x_els[i];	//!!!! To skip hillslope !!!!
            //for(j=0;j<assim_dim;j++)
            //	x_start[above_gauges[i]*assim_dim+j] += x_els[i*assim_dim+j];
        }
        VecRestoreArray(ws->x, &buffer);
    }

    //Send solution to everyone
    MPI_Bcast(x_start, allstates, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (verbose && my_rank == 0)
    {
        //unsigned int idxm[num_obs*steps_to_use];
        //double temp_matptr[(num_obs*steps_to_use*allstates_needed > allstates_needed*allstates_needed) ? num_obs*steps_to_use*allstates_needed : allstates_needed*allstates_needed];
        //for(i=0;i<num_obs*steps_to_use;i++)
        //    idxm[i] = i;

        //printf("x_start\n");
        //for(i=0;i<allstates;i++)
        //    printf("%.15e ",x_start[i]);
        //printf("\n");

        double* buffer;
        printf("difference (x)\n");
        VecGetArray(ws->x, &buffer);
        for (i = 0; i < allstates_needed; i++)
            printf("%.2e ", buffer[i]);
        printf("\n");
        VecRestoreArray(ws->x, &buffer);

        printf("d\n");
        VecGetArray(d, &buffer);
        for (i = 0; i < num_total_obs; i++)
            printf("%.2e ", buffer[i]);
        printf("\n");
        VecRestoreArray(d, &buffer);

        //printf("HM\n");
        //MatGetValues(*HM,num_obs*max_or_steps,idxm,allstates_needed,cols_allstates_needed,temp_matptr);
        //for(i=0;i<num_obs*max_or_steps;i++)
        //{
           // for(j=0;j<allstates_needed;j++)
              //  printf("%.15e ",temp_matptr[i*allstates_needed + j]);
           // printf(";\n");
        //}

        //printf("HTH\n");
        //MatGetValues(*HTH,allstates_needed,cols_allstates_needed,allstates_needed,cols_allstates_needed,temp_matptr);
        //for(i=0;i<allstates_needed;i++)
        //{
           // for(j=0;j<allstates_needed;j++)
              //  printf("%.15e ",temp_matptr[i*allstates_needed + j]);
           // printf(";\n");
        //}
    }

    //if (verbose)
    //{
    //    //Get q's produced from analysis (for testing)
    //    if (my_rank == 0)
    //    {
    //        //printf("q before\n");
    //        //Print_VECTOR(q,num_total_obs);
    //    }
    //    double first_diff = compute_diff(d_els, q, num_total_obs);

    //    ResetSysLS(sys, N, globals, t_b, x_start, problem_dim, globals->num_forcings, asynch->my_data);
    //    for (i = 0; i < N; i++)
    //        if (assignments[i] == my_rank || getting[i])
    //            //ReadInitData(globals->global_params,sys[i].params,NULL,0,sys[i].list->head->y_approx,globals->type,sys[i].diff_start,sys[i].no_ini_start,sys[i].user,NULL);
    //            custom_model->InitializeEqs(globals->global_params, sys[i].params, NULL, 0, sys[i].list->head->y_approx, globals->type, sys[i].diff_start, sys[i].no_ini_start, sys[i].user, NULL); //!!!! Should all states be reset? !!!!
    //    for (i = 0; i < asynch->globals->num_forcings; i++)
    //    {
    //        if (asynch->forcings[i]->flag == 3)	//!!!! I think .mon and binary files need this too !!!!
    //            Asynch_Set_Forcing_State(asynch, i, t_b, asynch->forcings[i]->first_file, asynch->forcings[i]->last_file);
    //    }
    //    for (i = 0; i < max_or_steps; i++)
    //    {
    //        globals->maxtime = t_b + (i)* inc;
    //        if (i)	Asynch_Advance(asynch, 0);

    //        for (j = 0; j < num_obs; j++)
    //        {
    //            owner = assignments[obs_locs[j]];
    //            my_link = (owner == my_rank);

    //            if (my_link)
    //                q[i*num_obs + j] = sys[obs_locs[j]]->list->tail->y_approx.ve[0];

    //            MPI_Bcast(&(q[i*num_obs + j]), 1, MPI_DOUBLE, owner, MPI_COMM_WORLD);
    //        }
    //    }

    //    if (my_rank == 0)
    //    {
    //        //printf("q after\n");
    //        //Print_VECTOR(q,num_total_obs);

    //        double second_diff = compute_diff(d_els, q, num_total_obs);
    //        printf("\nDifferences between q and data are %e %e\n\n", first_diff, second_diff);
    //    }
    //}



    //Clean up
    VecDestroy(&d);	//!!!! Blah !!!!

    MPI_Barrier(MPI_COMM_WORLD);
    double stop = MPI_Wtime();
    if (my_rank == 0)	printf("Total time for linear least squares fit: %.0f\n", stop - start);

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

    //return pow(result, 0.5);
    return result;
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
        for (j = 0; j < n; j++)	printf("%.2e ", A[i][j]);
        printf(";\n");
    }
    printf("\n");
}

void Print_VECTOR(double* v, unsigned int dim)
{
    unsigned int i;

    for (i = 0; i < dim; i++)
        printf("[%d]: %.2e\n", i, v[i]);
    printf(";\n");
}


////Read into memory the times and discharges stored in the DB.
//double*** DownloadGaugeReadings(unsigned int start_time, unsigned int stop_time, unsigned int** id_to_loc, unsigned int N, unsigned int* numlinks, unsigned int** ids, unsigned int** locs, unsigned int** numsteps)
//{
//    unsigned int i, j/*,k*/;
//    double ***data = NULL;
//    char query[1028];
//
//    if (my_rank == 0)
//    {
//        printf("Downloading gauge data...\n");
//        PGresult *res;
//
//        //Hard coding!
//        //unsigned int outlet = 434478;  //Turkey River above French Hollow
//        unsigned int outlet = 434514;	//Turkey River at Garber
//        //unsigned int outlet = 307864;  //Half Squaw Creek
//        //unsigned int outlet = 292254;	//Squaw Creek at Ames
//        //ConnData* conninfo = CreateConnData("dbname=model_ifc host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
//        ConnData* conninfo = CreateConnData("dbname=arch_usgs host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
//        ConnectPGDB(conninfo);
//
//        //Get link ids of gauges
//        //sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
//		//  SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, sensor_data AS C \
//		//  WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND B.link_id = A.link_id AND B.link_id != 311903 ORDER BY A.link_id;",outlet);
//        //sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
//		//	SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, sensor_data AS C \
//		//	WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND B.link_id = A.link_id AND B.link_id != 311903 AND B.link_id != 301218 AND B.link_id != 305680 ORDER BY A.link_id;",outlet);
//        sprintf(query, "WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
//			SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, master_usgs_gauges AS C \
//			WHERE B.link_id = A.link_id AND A.id = C.ifis_id AND A.link_id != 421097 AND A.link_id != 434582 ORDER BY link_id;", outlet);
//        res = PQexec(conninfo->conn, query);
//        CheckResError(res, "getting list of bridge sensor link ids");
//        *numlinks = PQntuples(res);
//
//        //Allocate space
//        *ids = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
//        *locs = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
//        *numsteps = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
//        data = (double***)malloc(*numlinks * sizeof(double**));
//
//        //Store the ids
//        for (i = 0; i < *numlinks; i++)
//            (*ids)[i] = atoi(PQgetvalue(res, i, 0));
//
//        //Download data
//        for (i = 0; i < *numlinks; i++)
//        {
//            //sprintf(query,"SELECT (extract('epoch' FROM date_trunc('minute', C.time2)) - %u)/60 AS unix_time, getdischarges((dist_bottom - measured_dist)*0.0328084::real,A.id)*0.0283168 AS q \
//			//	FROM env_pois_adv AS A, sensor_data AS C \
//			//	WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND A.link_id = %u AND to_timestamp(%u) <= C.time2 AND C.time2 <= to_timestamp(%u) \
//			//	ORDER BY unix_time;",start_time,(*ids)[i],start_time,stop_time);
//            sprintf(query, "SELECT (unix_timestamp-%u)/60 AS t,discharge*0.0283168 AS q FROM env_pois_adv AS A, master_usgs_gauges AS B \
//				WHERE A.id = B.ifis_id AND A.link_id = %u AND %u <= B.unix_timestamp AND B.unix_timestamp <= %u ORDER BY t;", start_time, (*ids)[i], start_time, stop_time);
//            res = PQexec(conninfo->conn, query);
//            CheckResError(res, "downloading bridge sensor data");
//            (*numsteps)[i] = PQntuples(res);
//            data[i] = (double**)malloc((*numsteps)[i] * sizeof(double*));
//
//            for (j = 0; j < (*numsteps)[i]; j++)
//            {
//                data[i][j] = (double*)malloc(2 * sizeof(double));
//                data[i][j][0] = atof(PQgetvalue(res, j, 0));
//                data[i][j][1] = atof(PQgetvalue(res, j, 1));
//            }
//        }
//
//        //Find locations from ids
//        for (i = 0; i < *numlinks; i++)
//            (*locs)[i] = find_link_by_idtoloc((*ids)[i], id_to_loc, N);
//
//        //Cleanup
//        PQclear(res);
//        DisconnectPGDB(conninfo);
//        ConnData_Free(conninfo);
//    }
//
//    /*
//        //Send data to all procs
//        //!!!! Would be better if only data for links assigned to this proc were available !!!!
//        //!!!! Only keep data for links here !!!!
//        MPI_Bcast(numlinks,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//        if(my_rank != 0)
//        {
//            *ids = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
//            *locs = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
//            *numsteps = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
//            data = (double***) malloc(*numlinks*sizeof(double**));
//        }
//        MPI_Bcast(*ids,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//        MPI_Bcast(*locs,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//        MPI_Bcast(*numsteps,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//
//        for(i=0;i<*numlinks;i++)
//        {
//            if(my_rank != 0)
//                data[i] = (double**) malloc((*numsteps)[i]*sizeof(double*));
//            for(j=0;j<(*numsteps)[i];j++)
//            {
//                if(my_rank != 0)	data[i][j] = (double*) malloc(2*sizeof(double));
//                MPI_Bcast(data[i][j],2,MPI_DOUBLE,0,MPI_COMM_WORLD);
//            }
//        }
//    */
//
//    MPI_Bcast(numlinks, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    if (my_rank != 0)
//    {
//        *numsteps = NULL;
//        *ids = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
//        *locs = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
//    }
//    MPI_Bcast(*ids, *numlinks, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(*locs, *numlinks, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//
//    return data;
//}



//int Output_Linkid(double t, VEC y_i, VEC global_params, VEC params, int state, void* user)
//{
//    return ((Link*)user)->ID;
//}
//
//void Set_Output_User_LinkID(asynchsolver* asynch)
//{
//    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
//    Link** sys = asynch->sys;
//
//    for (i = 0; i<my_N; i++)
//        sys[my_sys[i]].output_user = (void*)sys[my_sys[i]];
//}
//
//int Output_Timestamp(double t, VEC y_i, VEC global_params, VEC params, int state, void* user)
//{
//    CustomParams* forecastparams = (CustomParams*)user;
//    return (int)(round(t * 60.0 + forecastparams->offset - forecastparams->simulation_time_with_data) + 0.1);
//}


////Custom parameters for forecasting ***********************************************************
//void Init_Output_User_forecastparams(asynchsolver* asynch)
//{
//    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
//    Link** sys = asynch->sys;
//
//    for (i = 0; i < my_N; i++)
//        sys[my_sys[i]].output_user = malloc(sizeof(CustomParams));
//}
//
//void Free_Output_User_forecastparams(asynchsolver* asynch)
//{
//    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
//    Link** sys = asynch->sys;
//
//    for (i = 0; i < my_N; i++)
//    {
//        free(sys[my_sys[i]].output_user);
//        sys[my_sys[i]].output_user = NULL;
//    }
//}
//
//void Set_Output_User_forecastparams(asynchsolver* asynch, unsigned int offset, unsigned int time_with_data)
//{
//    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
//    Link** sys = asynch->sys;
//    CustomParams* forecastparams;
//
//    for (i = 0; i < my_N; i++)
//    {
//        forecastparams = (CustomParams*)sys[my_sys[i]].output_user;
//        forecastparams->ID = sys[my_sys[i]].ID;
//        forecastparams->offset = offset;
//        forecastparams->simulation_time_with_data = time_with_data;
//    }
//}
//
//void Init_Output_PeakflowUser_Offset(asynchsolver* asynch)
//{
//    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
//    Link** sys = asynch->sys;
//
//    for (i = 0; i < my_N; i++)
//        sys[my_sys[i]].peakoutput_user = malloc(sizeof(unsigned int));
//}
//
//void Free_Output_PeakflowUser_Offset(asynchsolver* asynch)
//{
//    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
//    Link** sys = asynch->sys;
//
//    for (i = 0; i < my_N; i++)
//    {
//        free(sys[my_sys[i]].peakoutput_user);
//        sys[my_sys[i]].peakoutput_user = NULL;
//    }
//}
//
//void Set_Output_PeakflowUser_Offset(asynchsolver* asynch, unsigned int offset)
//{
//    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
//    Link** sys = asynch->sys;
//
//    for (i = 0; i < my_N; i++)
//        *(unsigned int*)(sys[my_sys[i]].peakoutput_user) = offset;
//}
//

