#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if !defined(_MSC_VER)
#define ASYNCH_SLEEP sleep
#else
#include <windows.h>
#define ASYNCH_SLEEP Sleep
#endif

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

//#if defined(HAVE_PETSC)
#include <petsc.h>
//#endif

#include "assim_ls_methods.h"
#include "asynch_interface.h"
#include "assim_models.h"

void MM_mult(double** A,double** B,double** C,unsigned int m,unsigned int inner,unsigned int n);
void VECTOR_Copy(double* u,double* v,unsigned int dim);
int LinearLeastSquares(void* ptr);
void Print_MATRIX(double** A,unsigned int m,unsigned int n);
void Print_VECTOR(double* v,unsigned int dim);
double*** DownloadGaugeReadings(unsigned int start_time,unsigned int stop_time,unsigned int** id_to_loc,unsigned int N,unsigned int* numlinks,unsigned int** ids,unsigned int** locs,unsigned int** numsteps);
double compute_diff(double* d,double* q,unsigned int size);

typedef struct
{
	double *d,*HM_els,t_b,*x_start,*HM_buffer;
	Vec *RHS,*x,*invupareas,*B,*R;
	Mat *HM,*HTH,*HMTR;
	asynchsolver* asynch;
	KSP *ksp;
	unsigned int problem_dim,allstates,allstates_needed,inc,max_steps_to_use,steps_to_use;
	unsigned int *data_locs,numdata,*above_gauges,num_above,*cols_allstates_needed,assim_dim,*vareq_shift,*inv_vareq_shift;
	int *entries_d;
} AppCtxFull;


int print_out = 0;
int my_rank;
int np;

int main(int argc,char **argv)
{
	//MPI Stuff
	char help[] = "Here is a help message.\n";
	PetscInitialize(&argc,&argv,NULL,help);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	//Parse input
	if(argc < 2)
	{
		if(my_rank == 0)
		{
			printf("Command line parameter required:  A universal variable file (.gbl).\n");
			printf("\n");
		}
		PetscFinalize();
		return 1;
	}

	//Asynch solver init

	//Declare variables
	time_t start,stop,q_start,q_stop;
	unsigned int i,j,k,l,m,n;
	RKMethod** AllMethods;
	Link* current;
	char additional[16];	//For output filename
	//srand(time(NULL));	//!!!! Is this needed? !!!!

	//Init asynch object and the river network
	asynchsolver* asynch = Asynch_Init(MPI_COMM_WORLD,&argc,&argv);

	//Model 15
	//Asynch_Custom_Model(asynch,&SetParamSizes_Assim,&ConvertParams_Assim,&InitRoutines_Assim,&Precalculations_Assim,&ReadInitData_Assim);
	//Model 254
	//Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254,&Precalculations_Assim_254,&ReadInitData_Assim_254);
	//Model 254, q
	Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_q,&Precalculations_Assim_254,&ReadInitData_Assim_254_q);
	//Model 254, q and s_p
	//Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_qsp,&Precalculations_Assim_254,&ReadInitData_Assim_254_qsp);
	//Model 254, q and s_t
	//Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_qst,&Precalculations_Assim_254,&ReadInitData_Assim_254_qst);
	//Asynch_Custom_Partitioning(asynch,&Partition_METIS_ByEqs);

	//Data assimilation parameters
	unsigned int numdata,*ids,*numsteps,*data_locs;
	double forecast_window = 1440.0;
	double inc = 15.0;
	unsigned int max_steps_to_use = 10;
	unsigned int least_squares_iters = 1;
	//For model 254
	//unsigned int problem_dim = 7;	//!!!! Generalize this !!!!
	//unsigned int assim_dim = 4;	//!!!! Generalize this !!!!
	//For model 254 trim
	unsigned int problem_dim = 4;	//!!!! Generalize this !!!!
	unsigned int assim_dim = 4;	//!!!! Generalize this !!!!

	//Initialize system
	if(my_rank == 0)	printf("Reading global file...\n");
	Asynch_Parse_GBL(asynch,argv[1]);
	if(my_rank == 0)	printf("Loading network...\n");
	Asynch_Load_Network(asynch);

	//unsigned int data_locs[] = {2};	//Locations of links with data	!!!! Make sure these are locations and not IDs !!!!
	//unsigned int data_locs[] = {1,3,4};	//!!!! This should come from ReadSolution. Or should be calculated at least. !!!!
	//data_locs = (unsigned int*) malloc(3*sizeof(unsigned int));
	//data_locs[0] = 1; data_locs[1] = 3; data_locs[2] = 4;
	//double*** truesolution = ReadSolution("TempData/assim/yobservations.dat",asynch->id_to_loc,N,&numdata,&ids,&data_locs,&numsteps);
	//double*** truesolution = ReadSolution("TempData/assim/testobservations.dat",asynch->id_to_loc,asynch->N,&numdata,&ids,&data_locs,&numsteps);
	//double*** truesolution = ReadSolution("TempData/assim/254testobservations.dat",asynch->id_to_loc,asynch->N,&numdata,&ids,&data_locs,&numsteps);
	//double*** truesolution = DownloadGaugeReadings(1402790400,1405382400,asynch->id_to_loc,N,&numdata,&ids,&data_locs,&numsteps);	//June 15 2014 to July 15 2014
	//double*** truesolution = DownloadGaugeReadings(1403654400,1405382400,asynch->id_to_loc,asynch->N,&numdata,&ids,&data_locs,&numsteps);	//June 25 2014 to July 15 2014
	double*** truesolution = DownloadGaugeReadings(1430448300,1432953900,asynch->id_to_loc,asynch->N,&numdata,&ids,&data_locs,&numsteps);	//May 1 2015 to May 30 2015 (minus 135 mins)

	if(my_rank == 0)	printf("Partitioning network...\n");
	Find_Upstream_Links(asynch,problem_dim,0,inc,max_steps_to_use,data_locs,numdata);
	Asynch_Partition_Network(asynch);
	Clean_Upstream_Links(asynch);
	if(my_rank == 0)	printf("Loading parameters...\n");
	Asynch_Load_Network_Parameters(asynch,0);
	if(my_rank == 0)	printf("Reading dam and reservoir data...\n");
	Asynch_Load_Dams(asynch);
	if(my_rank == 0)	printf("Setting up numerical error data...\n");
	Asynch_Load_Numerical_Error_Data(asynch);
	if(my_rank == 0)	printf("Initializing model...\n");
	Asynch_Initialize_Model(asynch);
	Setup_Errors(asynch,problem_dim);
	if(my_rank == 0)	printf("Loading initial conditions...\n");
	Asynch_Load_Initial_Conditions(asynch);
	if(my_rank == 0)	printf("Loading forcings...\n");
	Asynch_Load_Forcings(asynch);
	if(my_rank == 0)	printf("Loading output data information...\n");
	Asynch_Load_Save_Lists(asynch);
	if(my_rank == 0)	printf("Finalizing network...\n");
	Asynch_Finalize_Network(asynch);
	if(my_rank == 0)	printf("Calculating initial step sizes...\n");
	Asynch_Calculate_Step_Sizes(asynch);

	//Prepare output files
	Asynch_Prepare_Temp_Files(asynch);
	Asynch_Write_Current_Step(asynch);
	Asynch_Prepare_Peakflow_Output(asynch);
	Asynch_Prepare_Output(asynch);

	//Pull data from asynch
	unsigned int my_N = asynch->my_N,N = asynch->N,*my_sys = asynch->my_sys,**id_to_loc = asynch->id_to_loc,num_forcings = asynch->GlobalVars->num_forcings;
	int *assignments = asynch->assignments;
	Link** sys = asynch->sys;
	short int *getting = asynch->getting;
	UnivVars *GlobalVars = asynch->GlobalVars;
	model* custom_model = asynch->custom_model;

	//Set print_time to t_0
	Asynch_Reset_Temp_Files(asynch,sys[my_sys[0]]->last_t);

	//Initialize choices
	double t_b = 0.0;
	double t_f = Asynch_Get_Total_Simulation_Time(asynch);
	unsigned int dim = assim_dim;
	unsigned int allstates = dim * N;
	double x_b[allstates];
	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank)
		{
			for(j=0;j<assim_dim;j++)
				x_b[i*assim_dim + j] = sys[i]->list->tail->y_approx.ve[j];	//!!!! Need to be able to specify which states are used !!!!
		}
		MPI_Bcast(&(x_b[i*assim_dim]),assim_dim,MPI_DOUBLE,assignments[i],MPI_COMM_WORLD);
	}
	double analysis[allstates];
	//double step = .1;	//!!!! Should be pulled from GlobalVars !!!!

	//Other initializations
	unsigned int numstepstotake,steps_to_use;
	unsigned int iterations = round((t_f - t_b) / inc);
	double* d = calloc(numdata,sizeof(double));
	//unsigned int total_steps = round((inc - t_b) / step);
	//unsigned int length_t = round((t_f - t_b) / step) + 1;

	//unsigned int start_idx = 0;
	double* d_full = calloc(numdata*max_steps_to_use,sizeof(double));
	double* x_start = (double*) malloc(allstates*sizeof(double));	//Values used to start asynch solver in tao solvers

	//Find locations unaffected by gauges
	//short int* bool_above_gauges;
	//unsigned int* above_gauges;
	//unsigned int num_above = GaugeDownstream(asynch,&above_gauges,&bool_above_gauges,data_locs,numdata);
	//unsigned int allstates_needed = num_above * 2;	//For q and s_p
	//unsigned int allstates_needed = num_above;	//For q
	unsigned int *vareq_shift,*inv_vareq_shift;

	//Call model specific data assimilation routines
	//For Model 254
	//unsigned int allstates_needed = Setup_Fitting_Data_Model254(asynch,data_locs,numdata);
	//For Model 254 trim, q
	Setup_Fitting_Data_Model254_q(asynch,data_locs,numdata);
	//For Model 254 trim, q and s_p
	//Setup_Fitting_Data_Model254_qsp(asynch,data_locs,numdata);
	//For Model 254 trim, q and s_t
	//Setup_Fitting_Data_Model254_qst(asynch,data_locs,numdata);
	unsigned int allstates_needed = BuildStateShift(asynch,allstates,data_locs,numdata,&vareq_shift,&inv_vareq_shift);


printf("allstates_needed: %u allstates: %u\n",allstates_needed,allstates);

	//For linear least squares
	Vec RHS,x;
	Mat HM,HTH,HMTR;
	KSP ksp;
	double *HM_buffer = (double*) malloc(allstates_needed*sizeof(double));	//!!!! Blah !!!!
	int *cols_allstates_needed = NULL;
	if(my_rank == 0)
	{
		VecCreateSeq(MPI_COMM_SELF,allstates_needed,&RHS);
		VecCreateSeq(MPI_COMM_SELF,allstates_needed,&x);
		//MatCreateSeqDense(MPI_COMM_SELF,numdata*max_steps_to_use,allstates,NULL,&HM);
		MatCreateSeqDense(MPI_COMM_SELF,allstates_needed,allstates_needed,NULL,&HTH);
		MatCreateSeqDense(MPI_COMM_SELF,allstates_needed,numdata*max_steps_to_use,NULL,&HMTR);
		//double* HM_els = (double*) malloc(allstates*sizeof(double));
		cols_allstates_needed = (int*) malloc(allstates_needed*sizeof(int));
		for(i=0;i<allstates_needed;i++)	cols_allstates_needed[i] = i;
		KSPCreate(MPI_COMM_SELF,&ksp);
		KSPSetOperators(ksp,HTH,HTH);
		KSPSetFromOptions(ksp);	//!!!! Do I need this? !!!!
		//MatAssemblyBegin(HTH,MAT_FINAL_ASSEMBLY);
		//MatAssemblyEnd(HTH,MAT_FINAL_ASSEMBLY);
	}
	int* entries_d = (int*) malloc(max_steps_to_use*numdata*sizeof(int));
	for(i=0;i<max_steps_to_use*numdata;i++)	entries_d[i] = i;

	//Transfer upstream areas to all procs
	char* my_links_needed = (char*) calloc(N,sizeof(char)),*links_needed = (char*) calloc(N,sizeof(char));
	double* invupareas = (double*) calloc(N,sizeof(double));
	upstream_data* updata;
	for(i=0;i<N;i++)
	{
		//if(assignments[i] == my_rank)	invupareas[i] = 1.0/(sys[i]->params.ve[GlobalVars->area_idx]*1e3);
		if(assignments[i] == my_rank)	invupareas[i] = 1.0;
		MPI_Bcast(&(invupareas[i]),1,MPI_DOUBLE,assignments[i],MPI_COMM_WORLD);
	}
	for(i=0;i<numdata;i++)	//Links needed for fitting
	{
		if(assignments[data_locs[i]] == my_rank)
		{
			current = sys[data_locs[i]];
			updata = (upstream_data*) (current->user);
			my_links_needed[data_locs[i]] = 1;
			for(j=0;j<current->numparents;j++)
			{
				for(k=0;k<updata->num_upstream[j];k++)
					my_links_needed[updata->upstream[j][k]] = 1;
			}
		}
	}
	MPI_Allreduce(my_links_needed,links_needed,N,MPI_CHAR,MPI_LOR,MPI_COMM_WORLD);
	free(my_links_needed);

	//Build weight matrices

	//!!!! Assuming only q is changing !!!!
	unsigned int curr_idx = 0;
	Vec B,R;
	if(my_rank == 0)
	{
		VecCreateSeq(MPI_COMM_SELF,allstates_needed,&B);
		for(i=0;i<N;i++)
		{
			if(links_needed[i])
				VecSetValue(B,curr_idx++,invupareas[i],INSERT_VALUES);
		}
		VecAssemblyBegin(B);
		VecAssemblyEnd(B);

		VecCreateSeq(MPI_COMM_SELF,numdata*max_steps_to_use,&R);
		for(i=0;i<numdata;i++)
		{
			for(j=0;j<max_steps_to_use;j++)
				VecSetValue(R,j*numdata + i,invupareas[data_locs[i]],INSERT_VALUES);
		}
		VecAssemblyBegin(R);
		VecAssemblyEnd(R);
	}
	free(links_needed);

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

		VecCreateSeq(MPI_COMM_SELF,numdata*max_steps_to_use,&R);
		for(i=0;i<numdata;i++)
		{
			for(j=0;j<max_steps_to_use;j++)
				VecSetValue(R,i*max_steps_to_use + j,invupareas[data_locs[i]],INSERT_VALUES);
		}
		VecAssemblyBegin(R);
		VecAssemblyEnd(R);
	}
	free(links_needed);
*/
/*
double* tempy;
VecGetArray(B,&tempy);
printf("B is (%u)\n",allstates_needed);
Print_VECTOR(tempy,allstates_needed);
VecRestoreArray(B,&tempy);

VecGetArray(R,&tempy);
printf("R is\n");
Print_VECTOR(tempy,max_steps_to_use*numdata);
VecRestoreArray(R,&tempy);
*/
	//Scale by area
	//printf("Scaling init discharges by area...\n");
	//AdjustDischarges_Scale(asynch,data_locs,d_full,numdata,x_b,allstates);

	//Prep PetSC
	if(my_rank == 0)	printf("\nPrepping PetSc...\n");
	AppCtxFull user;
	user.HM = &HM;
	user.HMTR = &HMTR;
	user.B = &B;
	user.R = &R;
	//user.HM_els = HM_els;
	user.HM_els = NULL;
	user.HM_buffer = HM_buffer;
	user.HTH = &HTH;
	user.RHS = &RHS;
	user.cols_allstates_needed = cols_allstates_needed;
	user.entries_d = entries_d;
	user.ksp = &ksp;
	//user.H = H;
	user.d = d_full;
	user.x_start = x_start;
	user.x = &x;
	user.asynch = asynch;
	user.problem_dim = problem_dim;
	user.assim_dim = assim_dim;
	user.allstates = allstates;
	user.allstates_needed = allstates_needed;
	user.vareq_shift = vareq_shift;
	user.inv_vareq_shift = inv_vareq_shift;
	user.inc = inc;
	user.max_steps_to_use = max_steps_to_use;
	user.steps_to_use = steps_to_use;
	user.data_locs = data_locs;
	user.numdata = numdata;
	user.t_b = t_b;
	//user.above_gauges = above_gauges;
	//user.num_above = num_above;

/*
	//Vec Lowerbds,Upperbds,Guess;
	PetscInt indices[allstates_needed];
	for(i=0;i<num_above;i++)
	{
		for(j=0;j<problem_dim;j++)
			indices[i*problem_dim + j] = above_gauges[i]*problem_dim + j;
	}
*/

	//Print out some information
	unsigned int my_eqs = 0,total_eqs;
	for(i=0;i<my_N;i++)	my_eqs += sys[my_sys[i]]->dim;
	MPI_Reduce(&my_eqs,&total_eqs,1,MPI_UNSIGNED,MPI_SUM,0,MPI_COMM_WORLD);

	printf("[%i]: Good to go with %u links (%u eqs).\n",my_rank,my_N,my_eqs);
	if(my_rank == 0)
	{
		sleep(1);
		printf("\nNetwork has a total of %u links and %u equations.\n\n",N,total_eqs);
		printf("Making calculations...\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(1);
	start = time(NULL);

	for(k=0;k<iterations;k++)
	{

		if(my_rank == 0)
		{
			printf("\n\n*************************************\n");
			printf("Iteration %u/%u. Background = %e.\n",k,iterations,t_b);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		//steps_to_use = (k+1 < max_steps_to_use+1) ? k+1 : max_steps_to_use;
		steps_to_use = (k < max_steps_to_use) ? k+1 : max_steps_to_use;

		//Set the forecast window
		Asynch_Set_Total_Simulation_Time(asynch,(t_f > t_b + forecast_window) ? t_b + forecast_window : t_f);

		//Get the new observations
MPI_Barrier(MPI_COMM_WORLD);
time(&q_start);
		FindAllDischarges(truesolution,0.0 + (k)*inc,numdata,numsteps,d);
MPI_Barrier(MPI_COMM_WORLD);
time(&q_stop);
if(my_rank == 0)
printf("Time to get new discharges: %.0f\n",difftime(q_stop,q_start));
		//printf("Out\n");

		if(k >= max_steps_to_use)
		{
			for(j=0;j<(max_steps_to_use-1)*(numdata);j++)
				d_full[j] = d_full[j+numdata];
			//for(j=0;j<(max_steps_to_use-1)*(numdata);j++)
			//	d_full[j] = d_full[j+1];
			for(j=0;j<numdata;j++)
				d_full[(max_steps_to_use-1)*(numdata) + j] = d[j];
		}
		else
		{
			for(j=0;j<numdata;j++)
				d_full[(steps_to_use-1)*numdata + j] = d[j];
		}


if(print_out && my_rank == 0)
//if(my_rank == 0)
{
printf("d_full\n");
Print_VECTOR(d_full,steps_to_use*numdata);
printf("\n");
}



if(k == 0)	//!!!! Move this outside the loop. Can I just change minimizer? !!!!
{
/*
	if(print_out && my_rank == 0)
	{
		for(j=0;j<np;j++)
		{
			if(my_rank == 0)
			{
				printf("*************\n");
				printf("x_b before scale\n");
				Print_VECTOR(x_b,allstates);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		sleep(2);
	}
*/
	AdjustDischarges_Scale(asynch,data_locs,d_full,numdata,x_b,allstates,assim_dim);
/*
	if(print_out && my_rank == 0)
	{
		for(j=0;j<np;j++)
		{
			if(my_rank == 0)
			{
				printf("x_b after scale\n");
				Print_VECTOR(x_b,allstates);
				printf("*************\n");
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		sleep(2);
	}
*/

/*
printf("254\n%u\n0.0\n\n",N);
for(i=0;i<N;i++)
{
	printf("%u %e %e %e %e\n",sys[i]->ID,x_b[4*i],x_b[4*i+1],x_b[4*i+2],x_b[4*i+3]);
}
sleep(60);
*/
}


		//Set bounds
		//LoadBounds_percentage(&Lowerbds,&Upperbds,x_b,allstates_needed,0.75);
		//TaoSetVariableBounds(tao,Lowerbds,Upperbds);

		//Set the initial guess
		user.steps_to_use = steps_to_use;
//		for(i=0;i<allstates_needed;i++)
//			minimizer[i] = x_b[indices[i]];	// !!!! Need a better name for minimizer... !!!!
		//VecSetValues(Guess,allstates_needed,indices,x_b,INSERT_VALUES);
		for(i=0;i<N;i++)	//Copy in states that won't be affected by the optimization solver
		{
			//if(!bool_above_gauges[i])
			{
				for(j=0;j<assim_dim;j++)
					x_start[i*assim_dim+j] = analysis[i*assim_dim+j] = x_b[i*assim_dim+j];
			}
		}

		//Set bounds
		//LoadBounds_percentage(&Lowerbds,&Upperbds,minimizer,allstates_needed,0.75);	//Used to be x_b, but now that we're using allstates_needed...
		//TaoSetVariableBounds(tao,Lowerbds,Upperbds);

		//Setup HM matrix
		if(my_rank == 0 && k < max_steps_to_use + 1)	//!!!! Is this right? Could be wasteful... !!!!
		{
			MatCreateSeqDense(MPI_COMM_SELF,numdata*steps_to_use,allstates_needed,NULL,&HM);
			MatAssemblyBegin(HM,MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(HM,MAT_FINAL_ASSEMBLY);
		}


//for(i=0;i<allstates;i++)	//Just to stop valgrind from whining...
//x_start[i] = x_b[i];

if(print_out && my_rank == 0)
{
//printf("Going in, minimizer...\n");
//Print_VECTOR(minimizer,allstates_needed);
printf("Going in, x_start...\n");
Print_VECTOR(x_start,allstates);

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

		//Calculate analysis
		int success = 1;
//while(success)
{
		//if(k)
		if(k >= max_steps_to_use)
		{
for(i=0;i<least_squares_iters;i++)
{
MPI_Barrier(MPI_COMM_WORLD);
			success = LinearLeastSquares(&user);
MPI_Barrier(MPI_COMM_WORLD);
			//bound_analysis(x_b,x_start,allstates);
}
		}
		else	success = 0;

		//Copy x_start to analysis  !!!! This shouldn't happen. Try a pointer dance. Or maybe something better... !!!!
		//if(success)
		{
			for(i=0;i<allstates;i++)
				analysis[i] = x_start[i];
		}
/*
		else
		{
			printf("!!!! Not applying least squares fit. !!!!\n");
			for(i=0;i<allstates;i++)
				analysis[i] = x_b[i];
		}
*/
}

/*
		for(i=0;i<num_above;i++)
		{
			for(j=0;j<assim_dim;j++)
				analysis[above_gauges[i]*assim_dim+j] = x_start[i*assim_dim+j];
		}
*/

		//Clean up HM matrix
		if(my_rank == 0 && k < max_steps_to_use)
		{
			MatDestroy(&HM);
		}

/*
		CHKERRQ(TaoSolve(tao));
		TaoConvergedReason reason;
		PetscInt its;
		PetscReal f,gnorm,cnorm,xdiff;
		//TaoGetTerminationReason(tao,&reason);
		TaoGetSolutionStatus(tao,&its,&f,&gnorm,&cnorm,&xdiff,&reason);
		TaoView(tao,PETSC_VIEWER_STDOUT_SELF);
		printf("Values are its = %i f = %e gnrom = %e cnorm = %e xdiff = %e reason = %i\n",its,f,gnorm,cnorm,xdiff,reason);

		if(reason != -6)
		{
			for(i=0;i<num_above;i++)
			{
				for(j=0;j<assim_dim;j++)
					analysis[above_gauges[i]*assim_dim+j] = minimizer[i*assim_dim+j];
			}
		}
		else
		{
			for(i=0;i<allstates;i++)	analysis[i] = x_b[i];
		}
*/

if(print_out && my_rank == 0)
{
//if(k == 10)
{
printf("x_b\n");
Print_VECTOR(x_b,allstates);
printf("\n");

printf("analysis t_b = %e k = %u last assim time = %e\n",t_b,k,t_b + 5.0*(steps_to_use-1));
Print_VECTOR(analysis,allstates);
//getchar();
}
}

		//Switch to model without variational eqs
		//For model 254
		//Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Model_254,&Precalculations_Assim_254,&ReadInitData_Assim_254);
		//For model 254 trim
		Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Model_252,&Precalculations_Assim_254,&ReadInitData_Assim_254_q);	//!!!! I think this is ok... !!!!
		Asynch_Initialize_Model(asynch);

MPI_Barrier(MPI_COMM_WORLD);
time(&q_start);

		//Advance to get new background
		ResetSysLS(sys,N,GlobalVars,t_b,analysis,assim_dim,num_forcings,asynch->my_data);
		for(i=0;i<N;i++)				//!!!! Put this into ResetSysLS? !!!!
		{
			if(assignments[i] == my_rank || getting[i])
			{
				current = sys[i];
				custom_model->InitializeEqs(GlobalVars->global_params,current->params,NULL,0,current->list->head->y_approx,GlobalVars->type,current->diff_start,current->no_ini_start,current->user,NULL);	//!!!! Blah. Shouldn't have to solve these... !!!!
				//ReadInitData(GlobalVars->global_params,current->params,NULL,0,current->list->head->y_approx,GlobalVars->type,current->diff_start,current->no_ini_start,current->user,NULL);	//!!!! Blah. Shouldn't have to solve these... !!!!
			}
		}
		if(steps_to_use == max_steps_to_use)
		{
			GlobalVars->maxtime = t_b + inc;
			Asynch_Advance(asynch,1);

			//Extract background !!!! Isn't there a routine for this? (there needs to be...) !!!!
			for(j=0;j<N;j++)
			{
				if(assignments[j] == my_rank)
				{
					for(l=0;l<dim;l++)
						x_b[j*dim+l] = sys[j]->list->tail->y_approx.ve[l];
				}
				MPI_Bcast(&(x_b[j*dim]),dim,MPI_DOUBLE,assignments[j],MPI_COMM_WORLD);
			}
		}
		else
			VECTOR_Copy(analysis,x_b,allstates);

		//Make forecasts
		//GlobalVars->maxtime = t_f;
//double blah = GlobalVars->maxtime;
		Asynch_Set_Total_Simulation_Time(asynch,(t_f > t_b + forecast_window) ? t_b + forecast_window : t_f);
//printf("Forecasting to %.0f, t_b = %.0f (well, really %.0f. Previous stop was %.0f)\n",(t_f > t_b + forecast_window) ? t_b + forecast_window : t_f,t_b,GlobalVars->maxtime,blah);

		Asynch_Advance(asynch,1);
MPI_Barrier(MPI_COMM_WORLD);
time(&q_stop);
if(my_rank == 0)	printf("Time for forecast: %.0f\n",difftime(q_stop,q_start));

		//Switch back to model with variational eqs
		//For model 254
		//Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254,&Precalculations_Assim_254,&ReadInitData_Assim_254);
		//For model 254 trim, q
		Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_q,&Precalculations_Assim_254,&ReadInitData_Assim_254_q);
		//Model 254, q and s_p
		//Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_qsp,&Precalculations_Assim_254,&ReadInitData_Assim_254_qsp);
		//Model 254, q and s_t
		//Asynch_Custom_Model(asynch,&SetParamSizes_Assim_254,&ConvertParams_Assim_254,&InitRoutines_Assim_254_qst,&Precalculations_Assim_254,&ReadInitData_Assim_254_qst);
		Asynch_Initialize_Model(asynch);

		//Go to next time
		if(k >= max_steps_to_use-1)
		{
			//start_idx += total_steps;
			t_b += inc;
			user.t_b = t_b;
		}
		//else
		//{
			//start_idx = 0;			// !!!! Or do nothing at all? !!!!
		//}

		//Reset the temporary files
		if(my_rank == 0)	printf("Creating output file...\n");
MPI_Barrier(MPI_COMM_WORLD);
time(&q_start);
		sprintf(additional,"%u",k);	//Add the iteration number to the output files
		Asynch_Create_Output(asynch,additional);
		//Asynch_Set_Temp_Files(asynch,t_b,(void*)&t_b,0);
		Asynch_Reset_Temp_Files(asynch,t_b);
MPI_Barrier(MPI_COMM_WORLD);
time(&q_stop);
if(my_rank == 0)	printf("Time to create output: %.0f\n",difftime(q_stop,q_start));
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);

	if(my_rank == 0)	printf("\nTime for calculations: %f. Writing results to disk...\nCalculations complete! All done!\n",difftime(stop,start));

	//Clean up
	//free(bool_above_gauges);
	//free(above_gauges);
	free(d);
	free(d_full);
	free(x_start);
	free(vareq_shift);
	free(inv_vareq_shift);

	if(my_rank == 0)
	{
		MatDestroy(&HM);
		MatDestroy(&HTH);
		VecDestroy(&RHS);
		VecDestroy(&x);
		VecDestroy(&R);
		VecDestroy(&B);
		MatDestroy(&HMTR);
	}
	//free(HM_els);
	free(HM_buffer);
	free(cols_allstates_needed);
	free(entries_d);
	if(my_rank == 0)
		KSPDestroy(&ksp);

	free(invupareas);
	free(ids);
	free(data_locs);
	if(truesolution)
	{
		for(i=0;i<numdata;i++)
		{
			for(j=0;j<numsteps[i];j++)
				free(truesolution[i][j]);
			free(truesolution[i]);
		}
		free(truesolution);
	}
	if(numsteps)	free(numsteps);

	//Petsc clean up
	PetscFinalize();

	//Asynch clean up
	Free_Upstream_Links(asynch);
	Asynch_Delete_Temporary_Files(asynch);
	Asynch_Free(asynch);

	return 0;
}


//This computes the least squares fit assuming the background and analysis difference is linear in the innovations. 
//HM is (numdata*steps_to_use) X allstates_needed
//HM_els is 1 X allstates (i.e. a 1D array)
//HM_buffer is 1 X allstates_needed
int LinearLeastSquares(void* ptr)
{
	unsigned int i,j,k,m,n,l,counter;
	Link* current;
	double factor;
	time_t start,stop,start2;
	short int my_link;
	int owner;
	upstream_data* updata;

	//Unpack ptr
	AppCtxFull* user = (AppCtxFull*) ptr;
	asynchsolver* asynch = user->asynch;
	unsigned int N = asynch->N;
	Link** sys = asynch->sys;
	UnivVars* GlobalVars = asynch->GlobalVars;
	int* assignments = asynch->assignments;
	short int* getting = asynch->getting;
	Mat *HM = user->HM,*HTH = user->HTH,*HMTR = user->HMTR;
	Vec *RHS = user->RHS,d,*x = user->x,*B = user->B,*R = user->R;
	KSP *ksp = user->ksp;
	double *d_els = user->d;
	unsigned int *data_locs = user->data_locs,assim_dim = user->assim_dim;
	unsigned int problem_dim = user->problem_dim,allstates = user->allstates,inc = user->inc,max_steps_to_use = user->max_steps_to_use,steps_to_use = user->steps_to_use,numdata = user->numdata;
	int *cols_allstates_needed = user->cols_allstates_needed,*entries_d = user->entries_d;
	double t_b = user->t_b,*x_els;
	unsigned int allstates_needed = user->allstates_needed;
	double *RHS_els,*x_start = user->x_start,*HM_buffer = user->HM_buffer;
	unsigned int max_or_steps = (steps_to_use == max_steps_to_use) ? max_steps_to_use : steps_to_use;	//!!!! Good job Scott. Good job... !!!!
	double q[max_or_steps*numdata];
	model* custom_model = asynch->custom_model;
	unsigned int *vareq_shift = user->vareq_shift,*inv_vareq_shift = user->inv_vareq_shift;

	time(&start);

	//Build a vector structure for d !!!! Obviously, this needs to not happen... !!!!
	//if(max_or_steps > allstates_needed)	printf("[%i]: Error: max_or_steps > allstates_needed (%u > %u)\n",my_rank,max_or_steps,allstates_needed);
	VecCreateSeq(MPI_COMM_SELF,max_or_steps*numdata,&d);	//!!!! This needs to be fixed !!!!
	VecSet(d,0.0);
	VecSetValues(d,max_or_steps*numdata,entries_d,d_els,INSERT_VALUES);
/*
//!!!! Try to scale the init conditions !!!!
AdjustDischarges_Scale(asynch,data_locs,d_els,numdata,x_start,allstates);
printf("Adjusted x_start to\n");
for(i=0;i<allstates;i++)	printf("%.15e ",x_start[i]);
printf("\n");
getchar();
*/
	//Initialize the system
	ResetSysLS(sys,N,GlobalVars,t_b,x_start,assim_dim,GlobalVars->num_forcings,asynch->my_data);
	for(i=0;i<N;i++)				//!!!! Put this into ResetSysLS? !!!!
		if(assignments[i] == my_rank || getting[i])
			custom_model->InitializeEqs(GlobalVars->global_params,sys[i]->params,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type,sys[i]->diff_start,sys[i]->no_ini_start,sys[i]->user,NULL); //!!!! Should all states be reset? !!!!
			//ReadInitData(GlobalVars->global_params,sys[i]->params,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type,sys[i]->diff_start,sys[i]->no_ini_start,sys[i]->user,NULL);	//!!!! Very inefficient. Too many checks. !!!!
	for(i=0;i<asynch->GlobalVars->num_forcings;i++)
	{
		if(asynch->forcings[i]->flag == 3)	//!!!! I think .mon and binary files need this too !!!!
			Asynch_Set_Forcing_State(asynch,i,t_b,asynch->forcings[i]->first_file,asynch->forcings[i]->last_file);
	}

	//Advance the system and extract the HM matrix
	for(i=0;i<max_or_steps;i++)	//!!!! Start at i=1? For i = 0, I don't think we need to set anything... !!!!		//HM here holds the values of M that are needed
	{
		GlobalVars->maxtime = t_b + (i) * inc;
		if(i)
		{
MPI_Barrier(MPI_COMM_WORLD);
time(&start2);
			Asynch_Advance(asynch,0);
MPI_Barrier(MPI_COMM_WORLD);
time(&stop);
if(my_rank == 0)
	printf("Time for advance %u: %.0f\n",i,difftime(stop,start2));
		}


//printf("ID = %u, t = %e\n",sys[data_locs[0]]->ID,sys[data_locs[0]]->last_t);
//Print_Vector(sys[data_locs[0]]->list->tail->y_approx);
//printf("**********\n");

		//Build HM
		for(j=0;j<numdata;j++)
		{
			current = sys[data_locs[j]];	//!!!! Assumes only discharges !!!!
			owner = assignments[data_locs[j]];
			my_link = (owner == my_rank);

			//From my link
			if(my_link)
			{
				updata = (upstream_data*) (current->user);
				//for(n=0;n<allstates;n++)	HM_els[n] = 0.0;

				//Pull out needed data
				for(n=0;n<allstates_needed;n++)
					HM_buffer[n] = 0.0;
				for(n=0;n<updata->num_fit_states;n++)
{
//printf("ID = %u | Loading %e (from %u) into spot %u\n",current->ID,current->list->tail->y_approx.ve[updata->fit_states[n]],updata->fit_states[n],vareq_shift[updata->fit_to_universal[n]]);
//sleep(1);
					HM_buffer[vareq_shift[updata->fit_to_universal[n]]] = current->list->tail->y_approx.ve[updata->fit_states[n]];
}

				//Extract calculationed q's (Just needed for testing. Maybe...)
				q[i*numdata+j] = current->list->tail->y_approx.ve[0];
			}

			MPI_Bcast(HM_buffer,allstates_needed,MPI_DOUBLE,owner,MPI_COMM_WORLD);	//!!!! Only proc 0 needs this !!!!
			MPI_Bcast(&(q[i*numdata+j]),1,MPI_DOUBLE,owner,MPI_COMM_WORLD);

			unsigned int row_idx = i*numdata + j;
//printf("Got %u, %u\n",allstates_needed,updata->num_fit_states);
//for(n=0;n<allstates_needed;n++)
//	printf("%e ",HM_buffer[n]);
//printf("\n");
//char r = getchar();
			if(my_rank == 0)
				MatSetValues(*HM,1,&row_idx,allstates_needed,cols_allstates_needed,HM_buffer,INSERT_VALUES);
		}
	}

/*
	//Zero out any unused rows
	if(i < max_steps_to_use)
	{
printf("Starting at %u, going to %u\n",i*numdata,max_steps_to_use*numdata);
		for(j=0;j<allstates;j++)
			HM_els[j] = 0.0;
		for(i=i*numdata;i<max_steps_to_use*numdata;i++)
			MatSetValues(*HM,1,&i,allstates,cols_allstates_needed,HM_els,INSERT_VALUES);
	}
*/

	//Assemble the HM matrix
	if(my_rank == 0)
	{
		MatAssemblyBegin(*HM,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(*HM,MAT_FINAL_ASSEMBLY);

//MPI_Barrier(MPI_COMM_WORLD);
time(&start2);
		//Calculate innovations
		double *d_tmp;
		VecGetArray(d,&d_tmp);
		for(i=0;i<max_or_steps*numdata;i++)	d_tmp[i] = d_tmp[i] - q[i];
		VecRestoreArray(d,&d_tmp);

		//Build the linear system
		//HMTR is allstates_needed x (numdata*max_or_steps)
		//HM is (numdata*max_or_steps) x allstates_needed
		MatTranspose(*HM,MAT_REUSE_MATRIX,HMTR);
		MatDiagonalScale(*HMTR,NULL,*R);
		MatMatMult(*HMTR,*HM,MAT_REUSE_MATRIX,PETSC_DEFAULT,HTH);
		MatDiagonalSet(*HTH,*B,ADD_VALUES);
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
		MatMult(*HMTR,d,*RHS);
		MatAssemblyBegin(*HTH,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(*HTH,MAT_FINAL_ASSEMBLY);

//MPI_Barrier(MPI_COMM_WORLD);
time(&stop);
if(my_rank == 0)
	printf("Time for matrix computations: %.0f\n",difftime(stop,start2));

	//Compute analysis
//MPI_Barrier(MPI_COMM_WORLD);
time(&start2);
		KSPSolve(*ksp,*RHS,*x);
//MPI_Barrier(MPI_COMM_WORLD);
time(&stop);
if(my_rank == 0)
	printf("Time for inversion: %.0f\n",difftime(stop,start2));

		//Copy new solution to x_start
		VecGetArray(*x,&x_els);
		//for(i=0;i<num_above;i++)
		for(i=0;i<allstates_needed;i++)	//!!!! I think this is right... !!!!
		{
//printf("i = %u %u\n",i,inv_vareq_shift[i]);
//sleep(1);
			x_start[inv_vareq_shift[i]] += x_els[i];
			//x_start[above_gauges[i]*assim_dim] += x_els[i];	//!!!! To skip hillslope !!!!
			//for(j=0;j<assim_dim;j++)
			//	x_start[above_gauges[i]*assim_dim+j] += x_els[i*assim_dim+j];
		}
		VecRestoreArray(*x,&x_els);
	}

	//Send solution to everyone
	MPI_Bcast(x_start,allstates,MPI_DOUBLE,0,MPI_COMM_WORLD);

if(print_out && my_rank == 0)
{
unsigned int idxm[numdata*max_steps_to_use];
double temp_matptr[(numdata*max_steps_to_use*allstates_needed > allstates_needed*allstates_needed) ? numdata*max_steps_to_use*allstates_needed : allstates_needed*allstates_needed];
for(i=0;i<numdata*max_steps_to_use;i++)	idxm[i] = i;

printf("x_start\n");
for(i=0;i<allstates;i++)	printf("%.15e ",x_start[i]);
printf("\n");

double* temp_ptr;
printf("difference (x)\n");
VecGetArray(*x,&temp_ptr);
for(i=0;i<allstates_needed;i++)	printf("%.15e ",temp_ptr[i]);
printf("\n");
VecRestoreArray(*x,&temp_ptr);

printf("d\n");
VecGetArray(d,&temp_ptr);
for(i=0;i<numdata*max_or_steps;i++)	printf("%.15e ",temp_ptr[i]);
printf("\n");

printf("HM\n");
MatGetValues(*HM,numdata*max_or_steps,idxm,allstates_needed,cols_allstates_needed,temp_matptr);
for(i=0;i<numdata*max_or_steps;i++)
{
	for(j=0;j<allstates_needed;j++)
		printf("%.15e ",temp_matptr[i*allstates_needed + j]);
	printf(";\n");
}

printf("HTH\n");
MatGetValues(*HTH,allstates_needed,cols_allstates_needed,allstates_needed,cols_allstates_needed,temp_matptr);
for(i=0;i<allstates_needed;i++)
{
	for(j=0;j<allstates_needed;j++)
		printf("%.15e ",temp_matptr[i*allstates_needed + j]);
	printf(";\n");
}

}

if(print_out)
{
//Get q's produced from analysis (for testing)
if(my_rank == 0)
{
//printf("q before\n");
//Print_VECTOR(q,max_or_steps*numdata);
}
double first_diff = compute_diff(d_els,q,max_or_steps*numdata);

ResetSysLS(sys,N,GlobalVars,t_b,x_start,problem_dim,GlobalVars->num_forcings,asynch->my_data);
for(i=0;i<N;i++)
	if(assignments[i] == my_rank || getting[i])
		//ReadInitData(GlobalVars->global_params,sys[i]->params,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type,sys[i]->diff_start,sys[i]->no_ini_start,sys[i]->user,NULL);
		custom_model->InitializeEqs(GlobalVars->global_params,sys[i]->params,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type,sys[i]->diff_start,sys[i]->no_ini_start,sys[i]->user,NULL); //!!!! Should all states be reset? !!!!
for(i=0;i<asynch->GlobalVars->num_forcings;i++)
{
	if(asynch->forcings[i]->flag == 3)	//!!!! I think .mon and binary files need this too !!!!
		Asynch_Set_Forcing_State(asynch,i,t_b,asynch->forcings[i]->first_file,asynch->forcings[i]->last_file);
}
for(i=0;i<max_or_steps;i++)
{
	GlobalVars->maxtime = t_b + (i) * inc;
	if(i)	Asynch_Advance(asynch,0);

	for(j=0;j<numdata;j++)
	{
		owner = assignments[data_locs[j]];
		my_link = (owner == my_rank);

		if(my_link)
			q[i*numdata+j] = sys[data_locs[j]]->list->tail->y_approx.ve[0];

		MPI_Bcast(&(q[i*numdata+j]),1,MPI_DOUBLE,owner,MPI_COMM_WORLD);
	}
}

if(my_rank == 0)
{
//printf("q after\n");
//Print_VECTOR(q,max_or_steps*numdata);

double second_diff = compute_diff(d_els,q,max_or_steps*numdata);
printf("\nDifferences between q and data are %e %e\n\n",first_diff,second_diff);
}
}



	//Clean up
	VecDestroy(&d);	//!!!! Blah !!!!

	MPI_Barrier(MPI_COMM_WORLD);
	time(&stop);
	if(my_rank == 0)	printf("Total time for linear least squares fit: %.0f\n",difftime(stop,start));

	return 0;
	//if(second_diff < first_diff)	return 1;
	//else				return 0;
}


double compute_diff(double* d,double* q,unsigned int size)
{
	int i;
	double result = 0.0;

	for(i=0;i<size;i++)
		result += (d[i] - q[i]) * (d[i] - q[i]);

	return pow(result,0.5);
}


//v = u
void VECTOR_Copy(double* u,double* v,unsigned int dim)
{
	unsigned int i;
	for(i=0;i<dim;i++)
		v[i] = u[i];
}

//C = A*B
//A= m x inner, B= inner x p, C= m x n
void MM_mult(double** A,double** B,double** C,unsigned int m,unsigned int inner,unsigned int n)
{
	unsigned int i,j,k;
	
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)	C[i][j] = 0.0;
		for(k=0;k<inner;k++)
		{
			for(j=0;j<n;j++)
				C[i][j] += A[i][k] * B[k][j];
		}
	}
}

void Print_MATRIX(double** A,unsigned int m,unsigned int n)
{
	unsigned int i,j;

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)	printf("%.15e ",A[i][j]);
		printf(";\n");
	}
	printf("\n");
}

void Print_VECTOR(double* v,unsigned int dim)
{
	unsigned int i;

	for(i=0;i<dim;i++)
		printf("%.15e ",v[i]);
	printf(";\n");
}


//Read into memory the times and discharges stored in a .dat file.
double*** DownloadGaugeReadings(unsigned int start_time,unsigned int stop_time,unsigned int** id_to_loc,unsigned int N,unsigned int* numlinks,unsigned int** ids,unsigned int** locs,unsigned int** numsteps)
{
	unsigned int i,j,k;
	double ***data = NULL;
	char query[1028];

	if(my_rank == 0)
	{
		printf("Downloading gauge data...\n");
		PGresult *res;

		//Hard coding!
		//unsigned int outlet = 434478;  //Turkey River above French Hollow
		unsigned int outlet = 434514;	//Turkey River at Garber
		//unsigned int outlet = 307864;  //Half Squaw Creek
		//unsigned int outlet = 292254;	//Squaw Creek at Ames
		ConnData* conninfo = CreateConnData("dbname=model_ifc host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
		//ConnData* conninfo = CreateConnData("dbname=arch_usgs host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
		ConnectPGDB(conninfo);

		//Get link ids of gauges
		//sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
			SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, sensor_data AS C \
			WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND B.link_id = A.link_id AND B.link_id != 311903 ORDER BY A.link_id;",outlet);
		//sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
			SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, sensor_data AS C \
			WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND B.link_id = A.link_id AND B.link_id != 311903 AND B.link_id != 301218 AND B.link_id != 305680 ORDER BY A.link_id;",outlet);
		//sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
			SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, master_usgs_gauges AS C \
			WHERE B.link_id = A.link_id AND A.id = C.ifis_id AND A.link_id != 421097 AND A.link_id != 434582 ORDER BY link_id;",outlet);
		sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
			SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, _data_usgs AS C \
			WHERE B.link_id = A.link_id AND A.id = C.ifis_id AND A.link_id != 418967 ORDER BY link_id;",outlet);	//Turkey River at Garber
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"getting list of bridge sensor link ids");
		*numlinks = PQntuples(res);

		//Allocate space
		*ids = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
		*locs = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
		*numsteps = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
		data = (double***) malloc(*numlinks*sizeof(double**));

		//Store the ids
		for(i=0;i<*numlinks;i++)
			(*ids)[i] = atoi(PQgetvalue(res,i,0));

		//Download data
		for(i=0;i<*numlinks;i++)
		{
			//sprintf(query,"SELECT (extract('epoch' FROM date_trunc('minute', C.time2)) - %u)/60 AS unix_time, getdischarges((dist_bottom - measured_dist)*0.0328084::real,A.id)*0.0283168 AS q \
				FROM env_pois_adv AS A, sensor_data AS C \
				WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND A.link_id = %u AND to_timestamp(%u) <= C.time2 AND C.time2 <= to_timestamp(%u) \
				ORDER BY unix_time;",start_time,(*ids)[i],start_time,stop_time);
			//sprintf(query,"SELECT (unix_timestamp-%u)/60 AS t,discharge*0.0283168 AS q FROM env_pois_adv AS A, master_usgs_gauges AS B \
				WHERE A.id = B.ifis_id AND A.link_id = %u AND %u <= B.unix_timestamp AND B.unix_timestamp <= %u ORDER BY t;",start_time,(*ids)[i],start_time,stop_time);
			sprintf(query,"SELECT (unix_time-%u)/60 AS t,discharge*0.0283168 AS q FROM env_pois_adv AS A, _data_usgs AS B \
				WHERE A.id = B.ifis_id AND A.link_id = %u AND %u <= B.unix_time AND B.unix_time <= %u ORDER BY t;",start_time,(*ids)[i],start_time,stop_time);
			res = PQexec(conninfo->conn,query);
			CheckResError(res,"downloading bridge sensor data");
			(*numsteps)[i] = PQntuples(res);
			data[i] = (double**) malloc((*numsteps)[i]*sizeof(double*));

			for(j=0;j<(*numsteps)[i];j++)
			{
				data[i][j] = (double*) malloc(2*sizeof(double));
				data[i][j][0] = atof(PQgetvalue(res,j,0));
				data[i][j][1] = atof(PQgetvalue(res,j,1));
			}
		}

		//Find locations from ids
		for(i=0;i<*numlinks;i++)
			(*locs)[i] = find_link_by_idtoloc((*ids)[i],id_to_loc,N);

		//Cleanup
		PQclear(res);
		DisconnectPGDB(conninfo);
		ConnData_Free(conninfo);
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

	MPI_Bcast(numlinks,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	if(my_rank != 0)
	{
		*numsteps = NULL;
		*ids = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
		*locs = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
	}
	MPI_Bcast(*ids,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(*locs,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	return data;
}


