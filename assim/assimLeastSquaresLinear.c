#include "rkmethods.h"
#include "system.h"
#include "comm.h"
#include "riversys.h"
#include "processdata.h"
#include <time.h>
#include <libpq-fe.h>
#include <mpi.h>
#include "misc.h"
#include "rainfall.h"
#include "assimlsmethods.h"
#include "solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <petsctao.h>
#include <petsc.h>
#include "asynch_interface.h"

void MM_mult(double** A,double** B,double** C,unsigned int m,unsigned int inner,unsigned int n);
void VECTOR_Copy(double* u,double* v,unsigned int dim);
void bound_analysis(double* x_b,double* analysis,unsigned int allstates);
int LinearLeastSquares(void* ptr);
int LinearLeastSquares2(void* ptr);
//PetscErrorCode EvaluateFunctionAndGradient(Tao tao, Vec X,PetscReal *f, Vec G, void *ptr);
//PetscErrorCode EvaluateHessian(Tao tao,Vec x,Mat *H,Mat *Hpre,MatStructure *flag,void *ptr);
void Print_MATRIX(double** A,unsigned int m,unsigned int n);
void Print_VECTOR(double* v,unsigned int dim);
unsigned int GaugeDownstream(asynchsolver* asynch,unsigned int** above_gauges,short int** bool_above_gauges,unsigned int* gauges,unsigned int numdata);
int AdjustDischarges_Scale(asynchsolver* asynch,unsigned int* data_locs,double* d,unsigned int numdata,double* x,unsigned int allstates);
double*** DownloadGaugeReadings(unsigned int start_time,unsigned int stop_time,unsigned int** id_to_loc,unsigned int N,unsigned int* numlinks,unsigned int** ids,unsigned int** locs,unsigned int** numsteps);
double compute_diff(double* d,double* q,unsigned int size);

//PetscErrorCode MyMonitor(Tao tao,void* ptr);
//PetscErrorCode DestroyMonitor(void** dunno);

int print_out = 0;
unsigned int eval_counter = 0;

typedef struct
{
	double **H,*d,*HM_els,t_b,*x_start,*HM_buffer;
	Vec *RHS,*x;
	Mat *HM,*HTH;
	asynchsolver* asynch;
	KSP *ksp;
	unsigned int problem_dim,allstates,allstates_needed,inc,max_steps_to_use,steps_to_use,*data_locs,numdata,*above_gauges,num_above,*cols_allstates_needed;
	int *entries_d;
} AppCtxFull;

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
	time_t start,stop;
	unsigned int i,j,k,l,m,n;
	RKMethod** AllMethods;
	Link* current;
	char additional[16];	//For output filename

	//Init asynch object and the river network
	asynchsolver* asynch = Asynch_Init(MPI_COMM_WORLD);
	Asynch_Parse_GBL(asynch,argv[1]);
	Asynch_Load_System(asynch);

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

	//Set print_time to t_0
/*
	for(i=0;i<my_N;i++)
	{
		current = sys[my_sys[i]];
		if(current->save_flag)
		{
			current->next_save = current->last_t;
			current->disk_iterations = 0;
		}
	}
*/
	Asynch_Reset_Temp_Files(asynch,sys[my_sys[0]]->last_t);
/*
	//Reserve space for backups, and initialize
	VEC** backup = (VEC**) malloc(N*sizeof(VEC*));
	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank || getting[i] == 1)
		{
			backup[i] = v_get(GlobalVars->dim);
			v_copy(sys[i]->list->tail->y_approx,backup[i]);
		}
		else
			backup[i] = NULL;
	}
*/
	//Least Squares init

	//Observations
	srand(time(NULL));
	unsigned int numdata,*ids,*numsteps,*data_locs;
	//unsigned int data_locs[] = {2};	//Locations of links with data	!!!! Make sure these are locations and not IDs !!!!
	//unsigned int data_locs[] = {1,3,4};	//!!!! This should come from ReadSolution. Or should be calculated at least. !!!!
	//data_locs = (unsigned int*) malloc(3*sizeof(unsigned int));
	//data_locs[0] = 1; data_locs[1] = 3; data_locs[2] = 4;
	//double*** truesolution = ReadSolution("TempData/assim/yobservations.dat",asynch->id_to_loc,N,&numdata,&ids,&data_locs,&numsteps);
	//double*** truesolution = ReadSolution("TempData/assim/testobservations.dat",asynch->id_to_loc,N,&numdata,&ids,&data_locs,&numsteps);
	//double*** truesolution = DownloadGaugeReadings(1402790400,1405382400,asynch->id_to_loc,N,&numdata,&ids,&data_locs,&numsteps);	//June 15 2014 to July 15 2014
	double*** truesolution = DownloadGaugeReadings(1403654400,1405382400,asynch->id_to_loc,N,&numdata,&ids,&data_locs,&numsteps);	//June 25 2014 to July 15 2014

	//Initialize choices
	unsigned int problem_dim = GlobalVars->problem_dim;
	double t_b = 0.0;
	double t_f = Asynch_Get_Total_Simulation_Time(asynch);
	double forecast_window = 1440.0;
	double inc = 15.0;
	unsigned int dim = problem_dim;
	unsigned int allstates = dim * N;
	double x_b[allstates];
	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank)
		{
			for(j=0;j<problem_dim;j++)
				x_b[i*problem_dim + j] = sys[i]->list->tail->y_approx->ve[j];
		}
		MPI_Bcast(&(x_b[i*problem_dim]),problem_dim,MPI_DOUBLE,assignments[i],MPI_COMM_WORLD);
	}
	double analysis[allstates];
	double step = .1;	//!!!! Should be pulled from GlobalVars !!!!
	double** H = (double**) malloc(numdata*sizeof(double*));
	for(i=0;i<numdata;i++)
		H[i] = (double*) calloc(allstates,sizeof(double));
	for(i=0;i<numdata;i++)
		H[i][data_locs[i]*dim] = 1.0;
	unsigned int max_steps_to_use = 10;

	//Other initializations
	unsigned int numstepstotake,steps_to_use;
	unsigned int iterations = round((t_f - t_b) / inc);
	double* d = calloc(numdata,sizeof(double));
	unsigned int total_steps = round((inc - t_b) / step);
	unsigned int length_t = round((t_f - t_b) / step) + 1;

	unsigned int start_idx = 0;
	double* d_full = calloc(numdata*max_steps_to_use,sizeof(double));
	double* x_start = (double*) malloc(allstates*sizeof(double));	//Values used to start asynch solver in tao solvers

	//Find locations unaffected by gauges
	short int* bool_above_gauges;
	unsigned int* above_gauges;
	unsigned int num_above = GaugeDownstream(asynch,&above_gauges,&bool_above_gauges,data_locs,numdata);
	//unsigned int allstates_needed = num_above * problem_dim;	//Number of states needed by the optimizer
	unsigned int allstates_needed = num_above;	//!!!! To skip hillslope !!!!
/*
	//!!!! For dropping low order links, and to skip hillslopes !!!!
	unsigned int allstates_needed = 0;
	short int* not_dropped = (short int*) calloc(N,sizeof(short int));
	for(i=0;i<numdata;i++)
	{
		current = sys[data_locs[i]];
//printf("Checking %u\n",current->ID);
		if(!not_dropped[data_locs[i]])
		{
			allstates_needed++;	//For current
			not_dropped[data_locs[i]] = 1;
			for(j=0;j<current->numparents;j++)
			{
				for(k=0;k<current->numupstream[j];k++)
				{
//printf("ID %u\n",sys[current->upstream[j][k]]->ID);
					if(!not_dropped[current->upstream[j][k]])
					{
						allstates_needed++;
						not_dropped[current->upstream[j][k]] = 1;
					}
					else
					{
//printf("breaking...\n");
						break;
					}
				}

				//allstates_needed += current->numupstream[j];	//For parents
				//for(k=0;k<current->numupstream[j];k++)	not_dropped[k] = 1;
			}
		}
	}
	unsigned int next_spot = 0;
	for(i=0;i<num_above;i++)
	{
		if(not_dropped[above_gauges[i]])	//Link was not dropped; keep in the list
		{
			above_gauges[next_spot++] = above_gauges[i];
		}
	}
	num_above = allstates_needed;
	above_gauges = (unsigned int*) realloc(above_gauges,num_above*sizeof(unsigned int));
	free(not_dropped);
*/

printf("allstates_needed: %u above gauges: %u\n",allstates_needed,num_above);
//for(i=0;i<num_above;i++)
//	printf("%u ",above_gauges[i]);
//printf("\n");

	//For linear least squares
	Vec RHS,x;
	Mat HM,HTH;
	VecCreateSeq(MPI_COMM_SELF,allstates_needed,&RHS);
	VecCreateSeq(MPI_COMM_SELF,allstates_needed,&x);
	//MatCreateSeqDense(MPI_COMM_SELF,numdata*max_steps_to_use,allstates,NULL,&HM);
	MatCreateSeqDense(MPI_COMM_SELF,allstates_needed,allstates_needed,NULL,&HTH);
	double* HM_els = (double*) malloc(allstates*sizeof(double));
	double *HM_buffer = (double*) malloc(allstates_needed*sizeof(double));	//!!!! Blah !!!!
	int* cols_allstates_needed = (int*) malloc(allstates_needed*sizeof(int));
	for(i=0;i<allstates_needed;i++)	cols_allstates_needed[i] = i;
	int* entries_d = (int*) malloc(max_steps_to_use*numdata*sizeof(int));
	for(i=0;i<max_steps_to_use*numdata;i++)	entries_d[i] = i;
	KSP ksp;
	KSPCreate(MPI_COMM_SELF,&ksp);
	KSPSetOperators(ksp,HTH,HTH);
	KSPSetFromOptions(ksp);	//!!!! Do I need this? !!!!
	//MatAssemblyBegin(HM,MAT_FINAL_ASSEMBLY);
	//MatAssemblyEnd(HM,MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(HTH,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(HTH,MAT_FINAL_ASSEMBLY);


	//Scale by area
	//printf("Scaling init discharges by area...\n");
	//AdjustDischarges_Scale(asynch,data_locs,d_full,numdata,x_b,allstates);

	//Prep PetSC
	if(my_rank == 0)	printf("\nPrepping PetSc...\n");
	AppCtxFull user;
	user.HM = &HM;
	user.HM_els = HM_els;
	user.HM_buffer = HM_buffer;
	user.HTH = &HTH;
	user.RHS = &RHS;
	user.cols_allstates_needed = cols_allstates_needed;
	user.entries_d = entries_d;
	user.ksp = &ksp;
	user.H = H;
	user.d = d_full;
	user.x_start = x_start;
	user.x = &x;
	user.asynch = asynch;
	user.problem_dim = problem_dim;
	user.allstates = allstates;
	user.allstates_needed = allstates_needed;
	user.inc = inc;
	user.max_steps_to_use = max_steps_to_use;
	user.steps_to_use = steps_to_use;
	user.data_locs = data_locs;
	user.numdata = numdata;
	user.t_b = t_b;
	user.above_gauges = above_gauges;
	user.num_above = num_above;

	//Vec Lowerbds,Upperbds,Guess;
	PetscInt indices[allstates_needed];
	//double* minimizer;
	for(i=0;i<num_above;i++)
	{
		for(j=0;j<problem_dim;j++)
			indices[i*problem_dim + j] = above_gauges[i]*problem_dim + j;
	}
/*
printf("Indices:\n");
for(i=0;i<allstates_needed;i++)
printf("%i ",indices[i]);
printf("\n");
*/

	//Compute data assim model
	if(my_rank == 0)
	{
		printf("\nGood to go with a total of %u links.\n",N);
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
		//printf("In\n");
		FindAllDischarges(truesolution,0.0 + (k)*inc,numdata,numsteps,d);
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


if(print_out)
//if(my_rank == 0)
{
printf("d_full\n");
Print_VECTOR(d_full,steps_to_use*numdata);
printf("\n");
}



if(k == 0)	//!!!! Move this outside the loop. Can I just change minimizer? !!!!
{
	if(print_out)
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
	AdjustDischarges_Scale(asynch,data_locs,d_full,numdata,x_b,allstates);
	if(print_out)
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
/*
printf("%u\n3\n\n",N);
for(i=0;i<N;i++)
{
	printf("%u %e %e\n",sys[i]->ID,x_b[2*i],x_b[2*i+1]);
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
				for(j=0;j<problem_dim;j++)
					x_start[i*problem_dim+j] = analysis[i*problem_dim+j] = x_b[i*problem_dim+j];
			}
		}

		//Set bounds
		//LoadBounds_percentage(&Lowerbds,&Upperbds,minimizer,allstates_needed,0.75);	//Used to be x_b, but now that we're using allstates_needed...
		//TaoSetVariableBounds(tao,Lowerbds,Upperbds);

		//Setup HM matrix
		if(k < max_steps_to_use + 1)	//!!!! Is this right? Could be a little wasteful... !!!!
		{
			MatCreateSeqDense(MPI_COMM_SELF,numdata*steps_to_use,allstates_needed,NULL,&HM);
			MatAssemblyBegin(HM,MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(HM,MAT_FINAL_ASSEMBLY);
		}


//for(i=0;i<allstates;i++)	//Just to stop valgrind from whining...
//x_start[i] = x_b[i];

if(print_out)
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
for(i=0;i<1;i++)
{
			//success = LinearLeastSquares(&user);
			success = LinearLeastSquares2(&user);
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
			for(j=0;j<problem_dim;j++)
				analysis[above_gauges[i]*problem_dim+j] = x_start[i*problem_dim+j];
		}
*/

		//Clean up HM matrix
		if(k < max_steps_to_use)
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
				for(j=0;j<problem_dim;j++)
					analysis[above_gauges[i]*problem_dim+j] = minimizer[i*problem_dim+j];
			}
		}
		else
		{
			for(i=0;i<allstates;i++)	analysis[i] = x_b[i];
		}
*/

if(print_out)
{
//if(k == 10)
{
printf("x_b\n");
Print_VECTOR(x_b,allstates);
printf("\n");

printf("analysis t_b = %e k = %u start_idx = %u last assim time = %e\n",t_b,k,start_idx,t_b + 5.0*(steps_to_use-1));
Print_VECTOR(analysis,allstates);
//getchar();
}
}

		//Advance to get new background
		ResetSysLS(sys,N,GlobalVars,t_b,analysis,problem_dim,num_forcings,asynch->my_data);
		for(i=0;i<N;i++)				//!!!! Put this into ResetSysLS? !!!!
		{
			if(assignments[i] == my_rank || getting[i])
			{
				current = sys[i];
				ReadInitData(GlobalVars->global_params,current->params,current->iparams,NULL,0,current->list->head->y_approx,GlobalVars->type,GlobalVars->diff_start,GlobalVars->no_ini_start,NULL);	//!!!! Blah. Shouldn't have to solve these... !!!!
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
						x_b[j*dim+l] = sys[j]->list->tail->y_approx->ve[l];
				}
				MPI_Bcast(&(x_b[j*dim]),dim,MPI_DOUBLE,assignments[j],MPI_COMM_WORLD);
			}
		}
		else
			VECTOR_Copy(analysis,x_b,allstates);

		//Make predictions
		//GlobalVars->maxtime = t_f;
double blah = GlobalVars->maxtime;
		Asynch_Set_Total_Simulation_Time(asynch,(t_f > t_b + forecast_window) ? t_b + forecast_window : t_f);
//printf("Forecasting to %.0f, t_b = %.0f (well, really %.0f. Previous stop was %.0f)\n",(t_f > t_b + forecast_window) ? t_b + forecast_window : t_f,t_b,GlobalVars->maxtime,blah);

time_t q_start,q_stop;
time(&q_start);
		Asynch_Advance(asynch,1);
MPI_Barrier(MPI_COMM_WORLD);
time(&q_stop);
if(my_rank == 0)	printf("Time for forecast: %.0f\n",difftime(q_stop,q_start));


		//Go to next time
		if(k >= max_steps_to_use-1)
		{
			start_idx += total_steps;
			t_b += inc;
			user.t_b = t_b;
		}
		else
		{
			start_idx = 0;			// !!!! Or do nothing at all? !!!!
		}

		//Reset the temporary files
		if(my_rank == 0)	printf("Creating output file...\n");
		sprintf(additional,"%u",k);	//Add the iteration number to the output files
		Asynch_Create_Output(asynch,additional);
		//Asynch_Set_Temp_Files(asynch,t_b,(void*)&t_b,0);
		Asynch_Reset_Temp_Files(asynch,t_b);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);

	if(my_rank == 0)	printf("\nTime for calculations: %f. Writing results to disk...\nCalculations complete! All done!\n",difftime(stop,start));

	//Clean up
	free(bool_above_gauges);
	free(above_gauges);
	free(d);
	free(d_full);
	for(i=0;i<numdata;i++)
		free(H[i]);
	free(H);
	free(x_start);

	MatDestroy(&HM);
	MatDestroy(&HTH);
	VecDestroy(&RHS);
	VecDestroy(&x);
	free(HM_els);
	free(HM_buffer);
	free(cols_allstates_needed);
	free(entries_d);
	KSPDestroy(&ksp);

	free(ids);
	free(data_locs);
	for(i=0;i<numdata;i++)
	{
		for(j=0;j<numsteps[i];j++)
			free(truesolution[i][j]);
		free(truesolution[i]);
	}
	free(truesolution);
	free(numsteps);

	//Petsc clean up
	//VecRestoreArray(Guess,&minimizer);
	//VecDestroy(&Lowerbds);
	//VecDestroy(&Upperbds);
	//VecDestroy(&Guess);
	//TaoDestroy(&tao);
	PetscFinalize();

	//Asynch clean up
	Asynch_Delete_Temporary_Files(asynch);
	Asynch_Free(asynch);

	return 0;
}


//This computes the least squares fit assuming the background and analysis difference is linear in the innovations. 
//HM is (numdata*steps_to_use) X allstates_needed
//HM_els is 1 X allstates (i.e. a 1D array)
//HM_buffer is 1 X allstates_needed
int LinearLeastSquares2(void* ptr)
{
	unsigned int i,j,k,m,n,l,counter;
	Link* current;
	double factor;
	time_t start,stop;
	short int my_link;
	int owner;

	//Unpack ptr
	AppCtxFull* user = (AppCtxFull*) ptr;
	asynchsolver* asynch = user->asynch;
	unsigned int N = asynch->N;
	Link** sys = asynch->sys;
	UnivVars* GlobalVars = asynch->GlobalVars;
	int* assignments = asynch->assignments;
	short int* getting = asynch->getting;
	Mat *HM = user->HM,*HTH = user->HTH;
	Vec *RHS = user->RHS,d,*x = user->x;
	KSP *ksp = user->ksp;
	double *d_els = user->d;
	unsigned int *above_gauges = user->above_gauges,num_above = user->num_above,*data_locs = user->data_locs;
	unsigned int problem_dim = user->problem_dim,allstates = user->allstates,inc = user->inc,max_steps_to_use = user->max_steps_to_use,steps_to_use = user->steps_to_use,numdata = user->numdata;
	int *cols_allstates_needed = user->cols_allstates_needed,*entries_d = user->entries_d;
	double t_b = user->t_b;
	unsigned int allstates_needed = user->allstates_needed;
	double *HM_els = user->HM_els,*RHS_els,*x_start = user->x_start,*HM_buffer = user->HM_buffer;
	unsigned int max_or_steps = (steps_to_use == max_steps_to_use) ? max_steps_to_use : steps_to_use;
	double q[max_or_steps*numdata];

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
*/

	//Initialize the system
	ResetSysLS(sys,N,GlobalVars,t_b,x_start,problem_dim,GlobalVars->num_forcings,asynch->my_data);
	for(i=0;i<N;i++)				//!!!! Put this into ResetSysLS? !!!!
		if(assignments[i] == my_rank || getting[i])
			ReadInitData(GlobalVars->global_params,sys[i]->params,sys[i]->iparams,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type,GlobalVars->diff_start,GlobalVars->no_ini_start,NULL);	//!!!! Very inefficient. Too many checks. !!!!
	for(i=0;i<asynch->GlobalVars->num_forcings;i++)
	{
		if(asynch->forcings[i]->flag == 3)	//!!!! I think .mon and binary files need this too !!!!
			Asynch_Set_Forcing_State(asynch,i,t_b,asynch->forcings[i]->first_file,asynch->forcings[i]->last_file);
	}

	//Advance the system and extract the HM matrix
	for(i=0;i<max_or_steps;i++)	//!!!! Start at i=1? For i = 0, I don't think we need to set anything... !!!!		//HM here holds the values of M that are needed
	{
		GlobalVars->maxtime = t_b + (i) * inc;
		if(i)	Asynch_Advance(asynch,0);

		//Build HM
		for(j=0;j<numdata;j++)
		{
			current = sys[data_locs[j]];	//!!!! Assumes only discharges !!!!
			owner = assignments[data_locs[j]];
			my_link = (owner == my_rank);

			//From my link
			if(my_link)
			{
				for(n=0;n<allstates;n++)	HM_els[n] = 0.0;
/*
				for(n=0;n<num_above;n++)	//!!!! Can this be used for better efficiency? !!!!
				{
					for(l=0;l<problem_dim;l++)
						HM_els[above_gauges[n]*problem_dim+l] = 0.0;
				}
*/

				for(n=0;n<problem_dim;n++)
					HM_els[current->location*problem_dim + n] = current->list->tail->y_approx->ve[problem_dim + (problem_dim-1) + n];

				//From parents
				counter = 3*problem_dim - 1;
				for(n=0;n<current->numparents;n++)
				{
					for(l=0;l<current->numupstream[n];l++)
					{
						for(m=0;m<problem_dim;m++)
							HM_els[current->upstream[n][l]*problem_dim + m] = current->list->tail->y_approx->ve[counter++];
					}
				}

				//Separate out the needed data
				for(n=0;n<num_above;n++)
				{
					HM_buffer[n] = HM_els[above_gauges[n]*problem_dim];	//!!!! To skip hillslope !!!!
					//for(l=0;l<problem_dim;l++)
					//	HM_buffer[n*problem_dim + l] = HM_els[above_gauges[n]*problem_dim+l];
				}


				//Extract calculationed q's (Just needed for testing. Maybe...)
				q[i*numdata+j] = current->list->tail->y_approx->ve[0];
			}

			MPI_Bcast(HM_buffer,allstates_needed,MPI_DOUBLE,owner,MPI_COMM_WORLD);
			MPI_Bcast(&(q[i*numdata+j]),1,MPI_DOUBLE,owner,MPI_COMM_WORLD);

			unsigned int row_idx = i*numdata + j;
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
	MatAssemblyBegin(*HM,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*HM,MAT_FINAL_ASSEMBLY);

	//Calculate innovations
	double *d_tmp;
	VecGetArray(d,&d_tmp);
	for(i=0;i<max_or_steps*numdata;i++)	d_tmp[i] = d_tmp[i] - q[i];
	VecRestoreArray(d,&d_tmp);

	//Build the linear system
	//HM is (numdata*max_or_steps) x allstates_needed
	MatTransposeMatMult(*HM,*HM,MAT_REUSE_MATRIX,PETSC_DEFAULT,HTH);
	for(i=0;i<allstates_needed;i++)
	{
		MatSetValue(*HTH,i,i,1.0,ADD_VALUES);  //!!!! To skip hillslope !!!!
		//if(i%2)	MatSetValue(*HTH,i,i,1e4,ADD_VALUES);
		//else	MatSetValue(*HTH,i,i,1.0,ADD_VALUES);
	}
	MatMultTranspose(*HM,d,*RHS);
	MatAssemblyBegin(*HTH,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*HTH,MAT_FINAL_ASSEMBLY);

	//Compute analysis
	KSPSolve(*ksp,*RHS,*x);

	//Copy new solution to x_start
	double* x_els;
	VecGetArray(*x,&x_els);
	for(i=0;i<num_above;i++)
	{
		x_start[above_gauges[i]*problem_dim] += x_els[i];	//!!!! To skip hillslope !!!!
		//for(j=0;j<problem_dim;j++)
		//	x_start[above_gauges[i]*problem_dim+j] += x_els[i*problem_dim+j];
	}
	VecRestoreArray(*x,&x_els);


if(print_out)
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
		ReadInitData(GlobalVars->global_params,sys[i]->params,sys[i]->iparams,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type,GlobalVars->diff_start,GlobalVars->no_ini_start,NULL);
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
			q[i*numdata+j] = sys[data_locs[j]]->list->tail->y_approx->ve[0];

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


//VecRestoreArray(d,&temp_ptr);

//}//End print_out

	//Clean up
	VecDestroy(&d);	//!!!! Blah !!!!

	MPI_Barrier(MPI_COMM_WORLD);
	time(&stop);
	if(my_rank == 0)	printf("Total time for linear least squares fit: %.0f\n",difftime(stop,start));

	return 0;
	//if(second_diff < first_diff)	return 1;
	//else				return 0;
}


//HM is (numdata*steps_to_use) X allstates_needed
//HM_els is 1 X allstates (i.e. a 1D array)
//HM_buffer is 1 X allstates_needed
int LinearLeastSquares(void* ptr)
{
	unsigned int i,j,k,m,n,l,counter;
	Link* current;
	double factor;
	time_t start,stop;

	//Unpack ptr
	AppCtxFull* user = (AppCtxFull*) ptr;
	asynchsolver* asynch = user->asynch;
	unsigned int N = asynch->N;
	Link** sys = asynch->sys;
	UnivVars* GlobalVars = asynch->GlobalVars;
	Mat *HM = user->HM,*HTH = user->HTH;
	Vec *RHS = user->RHS,d,*x = user->x;
	KSP *ksp = user->ksp;
	double *d_els = user->d;
	unsigned int *above_gauges = user->above_gauges,num_above = user->num_above,*data_locs = user->data_locs;
	unsigned int problem_dim = user->problem_dim,allstates = user->allstates,inc = user->inc,max_steps_to_use = user->max_steps_to_use,steps_to_use = user->steps_to_use,numdata = user->numdata;
	int *cols_allstates_needed = user->cols_allstates_needed,*entries_d = user->entries_d;
	double t_b = user->t_b;
	unsigned int allstates_needed = user->allstates_needed;
	double *HM_els = user->HM_els,*RHS_els,*x_start = user->x_start,*HM_buffer = user->HM_buffer;
	unsigned int max_or_steps = (steps_to_use == max_steps_to_use) ? max_steps_to_use : steps_to_use;
	double q[max_or_steps*numdata];

	time(&start);

	//Build a vector structure for d !!!! Obviously, this needs to not happen... !!!!
	//if(max_or_steps > allstates_needed)	printf("[%i]: Error: max_or_steps > allstates_needed (%u > %u)\n",my_rank,max_or_steps,allstates_needed);
printf("Setting d...\n");
	VecCreateSeq(MPI_COMM_SELF,max_or_steps*numdata,&d);	//!!!! This needs to be fixed !!!!
	VecSet(d,0.0);
	VecSetValues(d,max_or_steps*numdata,entries_d,d_els,INSERT_VALUES);
/*
//!!!! Try to scale the init conditions !!!!
AdjustDischarges_Scale(asynch,data_locs,d_els,numdata,x_start,allstates);
printf("Adjusted x_start to\n");
for(i=0;i<allstates;i++)	printf("%.15e ",x_start[i]);
printf("\n");
*/
	//Initialize the system
	ResetSysLS(sys,N,GlobalVars,t_b,x_start,problem_dim,GlobalVars->num_forcings,asynch->my_data);
	for(i=0;i<N;i++)				//!!!! Put this into ResetSysLS? !!!!
		ReadInitData(GlobalVars->global_params,sys[i]->params,sys[i]->iparams,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type,GlobalVars->diff_start,GlobalVars->no_ini_start,NULL);	//!!!! Very inefficient. Too many checks. !!!!
	for(i=0;i<asynch->GlobalVars->num_forcings;i++)
	{
		if(asynch->forcings[i]->flag == 3)	//!!!! I think .mon and binary files need this too !!!!
			Asynch_Set_Forcing_State(asynch,i,t_b,asynch->forcings[i]->first_file,asynch->forcings[i]->last_file);
	}

	//Advance the system and extract the HM matrix
	for(i=0;i<max_or_steps;i++)	//!!!! Start at i=1? For i = 0, I don't think we need to set anything... !!!!		//HM here holds the values of M that are needed
	{
		GlobalVars->maxtime = t_b + (i) * inc;
		if(i)	Asynch_Advance(asynch,0);

		//Build HM
		for(j=0;j<numdata;j++)
		{
			current = sys[data_locs[j]];	//!!!! Assumes only discharges !!!!
			for(n=0;n<allstates;n++)	HM_els[n] = 0.0;
			//for(n=0;n<allstates_needed;n++)	HM_els[n] = 0.0;
/*
			for(n=0;n<num_above;n++)	//!!!! Can this be used for better efficiency? !!!!
			{
				for(l=0;l<problem_dim;l++)
					HM_els[above_gauges[n]*problem_dim+l] = 0.0;
			}
*/

			//From my link
			for(n=0;n<problem_dim;n++)
				HM_els[current->location*problem_dim + n] = current->list->tail->y_approx->ve[problem_dim + (problem_dim-1) + n];

			//From parents
			counter = 3*problem_dim - 1;
			for(n=0;n<current->numparents;n++)
			{
				for(l=0;l<current->numupstream[n];l++)
				{
					for(m=0;m<problem_dim;m++)
						HM_els[current->upstream[n][l]*problem_dim + m] = current->list->tail->y_approx->ve[counter++];
				}
			}
/*
//Add extra weights
for(n=0;n<allstates;n+=2)
{
	//HM_els[n] *= 0.9;
	HM_els[n+1] *= 100.0;
}
*/

			//Separate out the needed data
			for(n=0;n<num_above;n++)
			{
				HM_buffer[n] = HM_els[above_gauges[n]*problem_dim];	//!!!! To skip hillslope !!!!
				//for(l=0;l<problem_dim;l++)
				//	HM_buffer[n*problem_dim + l] = HM_els[above_gauges[n]*problem_dim+l];
			}
			

//Extract calculationed q's (Just needed for testing. Maybe...)
q[i*numdata+j] = current->list->tail->y_approx->ve[0];

/*
printf("j = %u\n",j);
for(n=0;n<allstates;n++)
	printf("%.15e ",HM_els[n]);
printf(";\n");
for(n=0;n<allstates_needed;n++)
	printf("%.15e ",HM_buffer[n]);
printf(";\n");
*/
			unsigned int row_idx = i*numdata + j;
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
	MatAssemblyBegin(*HM,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*HM,MAT_FINAL_ASSEMBLY);


//!!!! New alteration to d !!!!
Vec v,x_0;
double *v_els,*d_tmp,*tempy;
VecCreateSeq(MPI_COMM_SELF,max_or_steps*numdata,&v);
VecGetArray(d,&d_tmp);
/*
printf("d_tmp\n");
for(i=0;i<max_or_steps*numdata;i++)	printf("%.15e ",d_tmp[i]);
printf("\n");
*/
VecCreateSeq(MPI_COMM_SELF,allstates_needed,&x_0);
VecGetArray(x_0,&tempy); //
for(n=0;n<num_above;n++) //
{
	tempy[n] = x_start[above_gauges[n]*problem_dim];	//!!!! To skip hillslope !!!!
	//for(l=0;l<problem_dim;l++)
	//	tempy[n*problem_dim + l] = x_start[above_gauges[n]*problem_dim+l];
}
VecRestoreArray(x_0,&tempy); //
//VecSetValues(x_0,allstates_needed,cols_allstates_needed,x_start,INSERT_VALUES);	//Old

/*
printf("x_0\n");
VecGetArray(x_0,&tempy);
for(i=0;i<allstates_needed;i++)	printf("%.15e ",tempy[i]);
printf("\n");
VecRestoreArray(x_0,&tempy);
*/

MatMult(*HM,x_0,v);
VecGetArray(v,&v_els);
/*
printf("v\n");
for(i=0;i<max_or_steps*numdata;i++)	printf("%.15e ",v_els[i]);
printf("\n");
*/
for(i=0;i<max_or_steps*numdata;i++)	d_tmp[i] = d_tmp[i] - q[i] + v_els[i];
VecRestoreArray(v,&v_els);
VecRestoreArray(d,&d_tmp);



/*
MatGetSize(*HM,&m,&n);
printf("HM is %u by %u\n",m,n);
MatGetSize(*HTH,&m,&n);
printf("HTH is %u by %u\n",m,n);
VecGetSize(*RHS,&m);
printf("RHS is %u\n",m);
VecGetSize(d,&m);
printf("d is %u\n",m);
VecGetSize(*x,&m);
printf("x is %u\n",m);
sleep(1);
*/
//printf("Here\n");
	//Build the linear system
	//HM is (numdata*max_or_steps) x allstates
	MatTransposeMatMult(*HM,*HM,MAT_REUSE_MATRIX,PETSC_DEFAULT,HTH);
	MatMultTranspose(*HM,d,*RHS);

//Modify the system
double HTH_vals[allstates_needed];
int row;
VecGetArray(*RHS,&tempy);
/*
//For yconn
row = 0;
HTH_vals[0] = 1.0;
HTH_vals[1] = -sys[0]->params->ve[2] / sys[1]->params->ve[2];
HTH_vals[2] = 0.0;
tempy[row] = 0.0;
*/

/*
//For test basin
row = 8;
HTH_vals[0] = 0.0;HTH_vals[1] = 0.0;HTH_vals[2] = 0.0;HTH_vals[3] = 0.0;HTH_vals[4] = 0.0;HTH_vals[5] = 0.0;HTH_vals[6] = 0.0;
HTH_vals[7] = 1.0;
HTH_vals[8] = -sys[9]->params->ve[2] / sys[10]->params->ve[2];
MatSetValues(*HTH,1,&row,allstates_needed,cols_allstates_needed,HTH_vals,INSERT_VALUES);
tempy[row] = 0.0;
/*
//Smooth out links immediately upstream from a gauge (unless that link has a gauge)
row = 1; //link 2
HTH_vals[0] = 0.0;HTH_vals[1] = 0.0;HTH_vals[2] = 0.0;HTH_vals[3] = 0.0;HTH_vals[4] = 0.0;HTH_vals[5] = 0.0;HTH_vals[6] = 0.0;HTH_vals[7] = 0.0;HTH_vals[8] = 0.0;
HTH_vals[0] = 1.0;
HTH_vals[1] = -10.0*sys[2]->params->ve[2] / sys[1]->params->ve[2];
MatSetValues(*HTH,1,&row,allstates_needed,cols_allstates_needed,HTH_vals,INSERT_VALUES);
tempy[row] = 0.0;

row = 6; //link 8
HTH_vals[0] = 0.0;HTH_vals[1] = 0.0;HTH_vals[2] = 0.0;HTH_vals[3] = 0.0;HTH_vals[4] = 0.0;HTH_vals[5] = 0.0;HTH_vals[6] = 0.0;HTH_vals[7] = 0.0;HTH_vals[8] = 0.0;
HTH_vals[0] = 1.0;
HTH_vals[6] = -10.0*sys[8]->params->ve[2] / sys[1]->params->ve[2];
MatSetValues(*HTH,1,&row,allstates_needed,cols_allstates_needed,HTH_vals,INSERT_VALUES);
tempy[row] = 0.0;

row = 5; //link 7
HTH_vals[0] = 0.0;HTH_vals[1] = 0.0;HTH_vals[2] = 0.0;HTH_vals[3] = 0.0;HTH_vals[4] = 0.0;HTH_vals[5] = 0.0;HTH_vals[6] = 0.0;HTH_vals[7] = 0.0;HTH_vals[8] = 0.0;
HTH_vals[3] = 1.0;
HTH_vals[5] = -10.0*sys[7]->params->ve[2] / sys[4]->params->ve[2];
MatSetValues(*HTH,1,&row,allstates_needed,cols_allstates_needed,HTH_vals,INSERT_VALUES);
tempy[row] = 0.0;
*/

MatAssemblyBegin(*HTH,MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(*HTH,MAT_FINAL_ASSEMBLY);
VecRestoreArray(*RHS,&tempy);

	//Compute analysis
	KSPSolve(*ksp,*RHS,*x);
//printf("Here\n");
	//Copy new solution to x_start
	double* x_els;
	VecGetArray(*x,&x_els);
	for(i=0;i<num_above;i++)
	{
		x_start[above_gauges[i]*problem_dim] = x_els[i];	//!!!! To skip hillslope !!!!
		//for(j=0;j<problem_dim;j++)
		//	x_start[above_gauges[i]*problem_dim+j] = x_els[i*problem_dim+j];
	}
	VecRestoreArray(*x,&x_els);

//printf("In...\n");
//getchar();

unsigned int idxm[numdata*max_steps_to_use];
double temp_matptr[(numdata*max_steps_to_use*allstates_needed > allstates_needed*allstates_needed) ? numdata*max_steps_to_use*allstates_needed : allstates_needed*allstates_needed];
for(i=0;i<numdata*max_steps_to_use;i++)	idxm[i] = i;

printf("x_start\n");
for(i=0;i<allstates;i++)	printf("%.15e ",x_start[i]);
printf("\n");

double* temp_ptr;
printf("x\n");
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


//Get q's produced from analysis (for testing)
printf("q before\n");
Print_VECTOR(q,max_or_steps*numdata);
double first_diff = compute_diff(d_els,q,max_or_steps*numdata);

ResetSysLS(sys,N,GlobalVars,t_b,x_start,problem_dim,GlobalVars->num_forcings,asynch->my_data);
for(i=0;i<N;i++)
	ReadInitData(GlobalVars->global_params,sys[i]->params,sys[i]->iparams,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type,GlobalVars->diff_start,GlobalVars->no_ini_start,NULL);
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
		q[i*numdata+j] = sys[data_locs[j]]->list->tail->y_approx->ve[0];
	}
}


printf("q after\n");
Print_VECTOR(q,max_or_steps*numdata);

double second_diff = compute_diff(d_els,q,max_or_steps*numdata);
printf("\nDifferences between q and data are %e %e\n\n",first_diff,second_diff);


VecRestoreArray(d,&temp_ptr);

//printf("Out\n");
//getchar();

	//Clean up
	VecDestroy(&d);	//!!!! Blah !!!!
	time(&stop);
	printf("Total time for linear least squares fit: %.0f\n",difftime(stop,start));


	//char r = getchar();
	//getchar();
	//return ('a' == r);

	if(second_diff < first_diff)	return 1;
	else				return 0;
}

double compute_diff(double* d,double* q,unsigned int size)
{
	int i;
	double result = 0.0;

	for(i=0;i<size;i++)
		result += (d[i] - q[i]) * (d[i] - q[i]);

	return pow(result,0.5);
}


void bound_analysis(double* x_b,double* analysis,unsigned int allstates)
{
	double lower_bd = 1e-6,upper_bd = -1.0,lower_bd_rel = 0.1,upper_bd_rel = 100.0;
	unsigned int i;

	//For discharges
	for(i=0;i<allstates;i+=2)
	{
		if(lower_bd_rel > 0.0 && analysis[i] < lower_bd_rel * x_b[i])	analysis[i] = lower_bd_rel * x_b[i];
		if(lower_bd > 0.0 && analysis[i] < lower_bd)		analysis[i] = lower_bd;
		if(upper_bd_rel > 0.0 && analysis[i] > upper_bd_rel * x_b[i])	analysis[i] = upper_bd_rel * x_b[i];
		if(upper_bd > 0.0 && analysis[i] > upper_bd)		analysis[i] = upper_bd;
	}
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

//Gives a list of all locations with a downstream gauge, sorted by location.
//This allocates an appropriate amount of space for above_gauges and bool_above_gauges.
//bool_above_gauges[i] == 1 if sys[i] is above a gauge. 0 if not.
//above_gauges[i] is the list of link locations above a gauge.
//Returns the total number of gauges above a gauge.
unsigned int GaugeDownstream(asynchsolver* asynch,unsigned int** above_gauges,short int** bool_above_gauges,unsigned int* gauges,unsigned int numdata)
{
	if(!above_gauges || !bool_above_gauges)	return 0;
	unsigned int N = asynch->N,num_above,i,j;
	(*bool_above_gauges) = (short int*) calloc(N,sizeof(short int));
	Link *current,*next,**sys = asynch->sys;

	//Set the location with gauges
	for(i=0;i<numdata;i++)
		(*bool_above_gauges)[gauges[i]] = 1;
	num_above = numdata;

	//Trickle down from each external link, until a gauge is reached, marking all links
	for(i=0;i<N;i++)
	{
		if(sys[i]->numparents == 0)
		{
			current = sys[i];
			next = current;

			//See if there is a gauge
			while(next && (*bool_above_gauges)[next->location] == 0)	next = next->c;

			if(next)	//Gauge found
			{
				for(;current != next;current = current->c)
				{
					(*bool_above_gauges)[current->location] = 1;
					num_above++;
				}
			}
		}
	}

	//Setup above_gauges
	*above_gauges = (unsigned int*) malloc(num_above*sizeof(unsigned int));

	j = 0;
	for(i=0;i<N;i++)
	{
		if((*bool_above_gauges)[i])
			(*above_gauges)[j++] = i;
	}

	return num_above;
}


//Scales discharge upstream by area from the gauges
//!!!! Going to need areas when calling this with multiple procs !!!!
int AdjustDischarges_Scale(asynchsolver* asynch,unsigned int* data_locs,double* d,unsigned int numdata,double* x,unsigned int allstates)
{
	//Unpack
	UnivVars* GlobalVars = asynch->GlobalVars;
	unsigned int N = asynch->N,**id_to_loc = asynch->id_to_loc;
	int *assignments = asynch->assignments;
	Link** sys = asynch->sys;
	unsigned int area_idx = GlobalVars->area_idx,problem_dim = GlobalVars->problem_dim;

	unsigned int i,j,loc,num_upstream,prev_loc,curr_loc;
	unsigned int *upstream = (unsigned int*) malloc(N*sizeof(unsigned int));
	short int *locs_set = (short int*) calloc(N,sizeof(short int));	//1 if the discharge is to be changed, 0 if not
	double *locs_newq = (double*) malloc(N*sizeof(double));	//The new value for the discharge at location i
	unsigned int *counter = (unsigned int*) malloc(N*sizeof(unsigned int));
	double* upstream_areas = (double*) malloc(N*sizeof(double));	//The upstream area of the link location i
	Link* current;
	double ratio;

	//Kinda crappy, but get the upstream area for each link
	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank)	upstream_areas[i] = sys[i]->params->ve[area_idx];
		MPI_Bcast(&(upstream_areas[i]),1,MPI_DOUBLE,assignments[i],MPI_COMM_WORLD);
	}

	//Initialize the new discharges to a negative value
	for(i=0;i<N;i++)	locs_newq[i] = -1.0;

	//Store gauge readings
	for(i=0;i<numdata;i++)
	{
		locs_set[data_locs[i]] = 1;
		locs_newq[data_locs[i]] = d[i];
	}

	//Set the discharges downstream from each gauge. This is for locations with no downstream gauages.
	//This follows each gauge downstream until a link is found with locs_set or an outlet.
	//counter is set to help track the closest gauge.
	for(i=0;i<N;i++)	counter[i] = N;

	for(i=0;i<numdata;i++)
	{
		loc = data_locs[i];
		current = sys[loc]->c;
		counter[loc] = 0;
		prev_loc = loc;

		if(current)
		{
			curr_loc = current->location;
			while(!locs_set[curr_loc])
			{
				if(counter[curr_loc] > counter[prev_loc]+1)
				{
					counter[curr_loc] = counter[prev_loc]+1;
					locs_newq[curr_loc] = upstream_areas[curr_loc] * locs_newq[loc] / upstream_areas[loc];
					//locs_newq[curr_loc] = current->params->ve[area_idx] * locs_newq[loc] / sys[loc]->params->ve[area_idx];
				}
				prev_loc = curr_loc;
				current = current->c;
				if(current)	curr_loc = current->location;
				else		break;
			}
		}
	}

	for(i=0;i<N;i++)
		if(counter[i] < N)	locs_set[i] = 1;

	//Perform the trickle down one last time. This is for links with no upstream or downstream gauges.
	for(i=0;i<N;i++)
	{
		if(!locs_set[i])
		{
			current = sys[i];
			num_upstream = 0;
			while(current && !(locs_set[current->location]))
			{
				upstream[num_upstream++] = current->location;
				current = current->c;
			}

			if(current)	//Gauge downstream
			{
				ratio = locs_newq[current->location] / upstream_areas[current->location];
				//ratio = locs_newq[current->location] / current->params->ve[area_idx];
				for(j=0;j<num_upstream;j++)
				{
					locs_set[upstream[j]] = 1;
					locs_newq[upstream[j]] = upstream_areas[upstream[j]] * ratio;
					//locs_newq[upstream[j]] = sys[upstream[j]]->params->ve[area_idx] * ratio;
				}
			}
		}
	}

	//Set the determined discharge. If a link's discharge was not determined by the above process, then it lies in a totally ungauged basin.
	for(i=0;i<N;i++)
	{
		if(locs_set[i])
			x[problem_dim*i] = locs_newq[i];
	}

	//Clean up
	free(locs_set);
	free(upstream);
	free(locs_newq);
	free(counter);
	free(upstream_areas);
	return 0;
}




//Read into memory the times and discharges stored in a .dat file.
double*** DownloadGaugeReadings(unsigned int start_time,unsigned int stop_time,unsigned int** id_to_loc,unsigned int N,unsigned int* numlinks,unsigned int** ids,unsigned int** locs,unsigned int** numsteps)
{
	unsigned int i,j,k;
	double*** data;
	char query[1028];

	if(my_rank == 0)
	{
		PGresult *res;

		//Hard coding!
		unsigned int outlet = 434514;	//Turkey River at Garber
		//unsigned int outlet = 307864;  //Half Squaw Creek
		//unsigned int outlet = 292254;	//Squaw Creek at Ames
		//ConnData* conninfo = CreateConnData("dbname=model_ifc host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
		ConnData* conninfo = CreateConnData("dbname=arch_usgs host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
		ConnectPGDB(conninfo);

		//Get link ids of gauges
		//sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
		//	SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, sensor_data AS C \
		//	WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND B.link_id = A.link_id AND B.link_id != 311903 ORDER BY A.link_id;",outlet);
		//sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
			SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, sensor_data AS C \
			WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND B.link_id = A.link_id AND B.link_id != 311903 AND B.link_id != 301218 AND B.link_id != 305680 ORDER BY A.link_id;",outlet);
		sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
			SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, master_usgs_gauges AS C \
			WHERE B.link_id = A.link_id AND A.id = C.ifis_id AND A.link_id != 421097 AND A.link_id != 434582 ORDER BY link_id;",outlet);
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
			sprintf(query,"SELECT unix_timestamp-%u AS t,discharge*0.0283168 AS q FROM env_pois_adv AS A, master_usgs_gauges AS B WHERE A.link_id = %u AND to_timestamp(%u) <= C.time2 AND C.time2 <= to_timestamp(%u) \
				ORDER BY t;",start_time,(*ids)[i],start_time,stop_time);
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


	//Send data to all procs
	//!!!! Would be better if only data for links assigned to this proc were available !!!!
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



	return data;
}

