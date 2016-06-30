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
#include <petsctao.h>
#include "asynch_interface.h"

void MM_mult(double** A,double** B,double** C,unsigned int m,unsigned int inner,unsigned int n);
void VECTOR_Copy(double* u,double* v,unsigned int dim);
void LoadBounds(Vec* Lowerbds,Vec* Upperbds,unsigned int size);
void LoadBounds_percentage(Vec* Lowerbds,Vec* Upperbds,double* state,unsigned int size,double percent);
PetscErrorCode EvaluateFunctionAndGradient(Tao tao, Vec X,PetscReal *f, Vec G, void *ptr);
//PetscErrorCode EvaluateHessian(Tao tao,Vec x,Mat *H,Mat *Hpre,MatStructure *flag,void *ptr);
void Print_MATRIX(double** A,unsigned int m,unsigned int n);
void Print_VECTOR(double* v,unsigned int dim);
unsigned int GaugeDownstream(asynchsolver* asynch,unsigned int** above_gauges,short int** bool_above_gauges,unsigned int* gauges,unsigned int numdata);
int AdjustDischarges_Scale(asynchsolver* asynch,unsigned int* data_locs,double* d,unsigned int numdata,double* x,unsigned int allstates);
double*** DownloadGaugeReadings(unsigned int start_time,unsigned int stop_time,unsigned int** id_to_loc,unsigned int N,unsigned int* numlinks,unsigned int** ids,unsigned int** locs,unsigned int** numsteps);

//PetscErrorCode MyMonitor(Tao tao,void* ptr);
//PetscErrorCode DestroyMonitor(void** dunno);

int print_out = 1;
unsigned int eval_counter = 0;

typedef struct
{
	double *HM_full,**H,*d,*x_start,t_b;
	asynchsolver* asynch;
	unsigned int problem_dim,allstates,allstates_needed,inc,max_steps_to_use,steps_to_use,*data_locs,numdata,*above_gauges,num_above,*cols_allstates;
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
	//unsigned int *my_sys,my_N,my_max_nodes,N,nummethods,my_save_size,save_size,i,j,k,l,m,n,*save_list,peaksave_size;
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
	unsigned int numlinks_obs,*ids,*numsteps,*data_locs;
	//unsigned int data_locs[] = {2};	//Locations of links with data	!!!! Make sure these are locations and not IDs !!!!
	//unsigned int data_locs[] = {1,3,4};	//!!!! This should come from ReadSolution. Or should be calculated at least. !!!!
	//data_locs = (unsigned int*) malloc(3*sizeof(unsigned int));
	//data_locs[0] = 1; data_locs[1] = 3; data_locs[2] = 4;
	//double*** truesolution = ReadSolution("TempData/assim/yobservations.dat",&numlinks_obs,&ids,&numsteps);
	//double*** truesolution = ReadSolution("TempData/assim/testobservations.dat",&numlinks_obs,&ids,&numsteps);
	double*** truesolution = ReadSolution("TempData/assim/testobservations.dat",asynch->id_to_loc,N,&numlinks_obs,&ids,&data_locs,&numsteps);
	//double*** truesolution = DownloadGaugeReadings(1402790400,1405382400,asynch->id_to_loc,N,&numlinks_obs,&ids,&data_locs,&numsteps);

	//Initialize choices
	unsigned int problem_dim = GlobalVars->problem_dim;
	double t_b = 0.0;
	double t_f = 300.0;
	double inc = 5.0;
	unsigned int numdata = numlinks_obs;	//!!!! Get rid of numlinks_obs !!!!
	unsigned int dim = problem_dim;
	unsigned int allstates = dim * N;
	//double x_b[6] = {10, 0, 10, 0, 10, 0};
	double x_b[allstates];
	for(i=0;i<N;i++)	for(j=0;j<problem_dim;j++)
		x_b[i*problem_dim + j] = sys[i]->list->tail->y_approx->ve[j];
	double analysis[allstates];
	//unsigned int data = 2*2+1;
	double step = .1;
	double** H = (double**) malloc(numdata*sizeof(double*));
	for(i=0;i<numdata;i++)
		H[i] = (double*) calloc(allstates,sizeof(double));
	//H[0][4] = 1.0;
	for(i=0;i<numdata;i++)
		H[i][data_locs[i]*dim] = 1.0;
	unsigned int max_steps_to_use = 4;

	//Other initializations
	unsigned int numstepstotake,steps_to_use;
	unsigned int iterations = round((t_f - t_b) / inc);
	double* d = calloc(numdata,sizeof(double));
	unsigned int total_steps = round((inc - t_b) / step);
	unsigned int length_t = round((t_f - t_b) / step) + 1;

	unsigned int start_idx = 0;
	double* d_full = calloc(numdata*max_steps_to_use,sizeof(double));
	double* HM_full = (double*) calloc(allstates,sizeof(double));	//!!!! For full least squares !!!!
	//double** HM_full = (double**) malloc(numdata*max_steps_to_use*sizeof(double*));
	//for(i=0;i<numdata*max_steps_to_use;i++)
	//	HM_full[i] = (double*) calloc(allstates,sizeof(double));
	double* x_start = (double*) malloc(allstates*sizeof(double));	//Values used to start asynch solver in tao solvers

	//Find locations unaffected by gauges
	short int* bool_above_gauges;
	unsigned int* above_gauges;
	unsigned int num_above = GaugeDownstream(asynch,&above_gauges,&bool_above_gauges,data_locs,numdata);
	unsigned int allstates_needed = num_above * problem_dim;	//Number of states needed by the optimizer
/*
printf("above gauges: %u\n",num_above);
for(i=0;i<num_above;i++)
	printf("%u ",above_gauges[i]);
printf("\n");
*/

	//Scale by area
	//printf("Scaling init discharges by area...\n");
	//AdjustDischarges_Scale(asynch,data_locs,d_full,numdata,x_b,allstates);

	//Prep TAO
	if(my_rank == 0)	printf("\nPrepping TAO...\n");
	AppCtxFull user;
	user.HM_full = HM_full;
	user.H = H;
	user.d = d_full;
	user.x_start = x_start;

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
	//user.workspace = asynch->workspace;

	Vec Lowerbds,Upperbds,Guess;
	PetscInt indices[allstates_needed];
	double* minimizer;
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
	Tao tao;
	TaoCreate(PETSC_COMM_SELF,&tao);
	TaoSetType(tao,TAOBLMVM);
	//TaoSetType(tao,"tao_blmvm");
	//TaoSetType(tao,"tao_tron");
	//TaoSetType(tao,"tao_gpcg");
	//TaoSetType(tao,"tao_bqpip");
	VecCreateSeq(MPI_COMM_SELF,allstates_needed,&Lowerbds);
	VecCreateSeq(MPI_COMM_SELF,allstates_needed,&Upperbds);
	VecCreateSeq(MPI_COMM_SELF,allstates_needed,&Guess);
	TaoSetInitialVector(tao,Guess);
	TaoSetObjectiveAndGradientRoutine(tao,EvaluateFunctionAndGradient,&user);
	VecGetArray(Guess,&minimizer);	//Grab solution

	//TaoSetTolerances(tao,1e-8,1e-8,1e-8,1e-8,0.0);
	//TaoSetTolerances(tao,1e-6,1e-6,1e-5,1e-5,0.0);
	TaoSetTolerances(tao,1e-3,1e-4,1e-3,1e-4,0.0);
	//TaoSetMaximumFunctionEvaluations(tao,100);
	TaoSetFromOptions(tao);

	//Compute data assim model
	if(my_rank == 0)
	{
		printf("\nGood to go with a total of %u links.\n",N);
		printf("Making calculations...\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	for(k=0;k<iterations;k++)
	{
		printf("\n\n*************************************\n");
		printf("Iteration %u/%u. Background = %e.\n",k,iterations,t_b);

		//steps_to_use = (k+1 < max_steps_to_use+1) ? k+1 : max_steps_to_use;
		steps_to_use = (k < max_steps_to_use) ? k+1 : max_steps_to_use;

		//Get the new observations
		printf("In\n");
		FindAllDischarges(truesolution,0.0 + (k)*inc,numlinks_obs,numsteps,d);
		printf("Out\n");
/*
		for(i=0;i<numdata;i++)
		{
			current = sys[data_locs[i]];
			//d[i] = FindDischarge(truesolution,current->ID,0.0 + (k+1)*inc,numlinks_obs,ids,numsteps,0.0);
			d[i] = FindDischarge(truesolution,current->ID,0.0 + (k)*inc,numlinks_obs,ids,numsteps,0.0);
		}
*/

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
{
//printf("Guess %u\n",steps_to_use);
//Print_VECTOR(x_b,allstates);
//printf("indices\n");
//for(i=0;i<allstates;i++)	printf("%i ",indices[i]);
//printf("\n");
//printf("HM\n");
//Print_MATRIX(HM_full,max_steps_to_use,allstates);
//printf("\n");
printf("d_full\n");
Print_VECTOR(d_full,steps_to_use*numdata);
printf("\n");
}

/*
if(k == 0)	//!!!! Move this outside the loop. Can I just change minimizer? !!!!
{
printf("*************\n");
printf("x_b before scale\n");
Print_VECTOR(x_b,allstates);
AdjustDischarges_Scale(asynch,data_locs,d_full,numdata,x_b,allstates);
printf("x_b after scale\n");
Print_VECTOR(x_b,allstates);
printf("*************\n");
}
*/
		//Set bounds
		LoadBounds_percentage(&Lowerbds,&Upperbds,x_b,allstates_needed,0.75);
		TaoSetVariableBounds(tao,Lowerbds,Upperbds);

		//Set the initial guess
		user.steps_to_use = steps_to_use;
		for(i=0;i<allstates_needed;i++)
			minimizer[i] = x_b[indices[i]];	// !!!! Need a better name for minimizer... !!!!
		//VecSetValues(Guess,allstates_needed,indices,x_b,INSERT_VALUES);
		for(i=0;i<N;i++)	//Copy in states that won't be affected by the optimization solver
		{
			if(!bool_above_gauges[i])
			{
				for(j=0;j<problem_dim;j++)
					x_start[i*problem_dim+j] = analysis[i*problem_dim+j] = x_b[i*problem_dim+j];
			}
		}

		//Set bounds
		LoadBounds_percentage(&Lowerbds,&Upperbds,minimizer,allstates_needed,0.75);	//Used to be x_b, but now that we're using allstates_needed...
		TaoSetVariableBounds(tao,Lowerbds,Upperbds);


for(i=0;i<allstates;i++)	//Just to stop valgrind from whining...
x_start[i] = x_b[i];


printf("Going in, minimizer...\n");
Print_VECTOR(minimizer,allstates_needed);
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


		//Calculate analysis
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

/*
			for(i=0;i<N;i++)
			{
				if(bool_above_gauges[i])
					for(j=0;j<problem_dim;j++)
						analysis[i*problem_dim+j] = minimizer[i*problem_dim+j];
				else
					for(j=0;j<problem_dim;j++)
						analysis[i*problem_dim+j] = x_b[i*problem_dim+j];
			}
*/
		}
		else
		{
			for(i=0;i<allstates;i++)	analysis[i] = x_b[i];
		}

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

printf("Here_ %i\n",k);

		//Advance to get new background
		ResetSysLS(sys,N,GlobalVars,t_b,analysis,problem_dim,num_forcings);
		for(i=0;i<N;i++)				//!!!! Put this into ResetSysLS? !!!!
			ReadInitData(GlobalVars->global_params,sys[i]->params,sys[i]->iparams,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type,GlobalVars->diff_start,GlobalVars->no_ini_start,NULL);	//!!!! Blah. Shouldn't have to solve these... !!!!
printf("Here %i\n",k);
		if(steps_to_use == max_steps_to_use)
		{
			GlobalVars->maxtime = t_b + inc;
			Asynch_Advance(asynch,1);

			//Extract background !!!! Isn't there a routine for this? !!!!
			for(j=0;j<N;j++)
				for(l=0;l<dim;l++)
					x_b[j*dim+l] = sys[j]->list->tail->y_approx->ve[l];
		}
		else
			VECTOR_Copy(analysis,x_b,allstates);
printf("Here %i\n",k);
		//Make predictions
		GlobalVars->maxtime = t_f;
		//AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,1,outputfile);

time_t q_start,q_stop;
time(&q_start);
		Asynch_Advance(asynch,1);
time(&q_stop);
printf("Time for forecast: %.0f\n",difftime(q_stop,q_start));
printf("Here %i\n",k);
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
printf("Here %i\n",k);
		//Reset the temporary files
		printf("Creating output file...\n");
		sprintf(additional,"%u",k);	//Add the iteration number to the output files
		//Process_Data(sys,GlobalVars,N,save_list,save_size,my_save_size,id_to_loc,assignments,NULL,additional);
		Asynch_Create_Output(asynch,additional);
		//SetTempFiles(t_b,sys,N,outputfile,GlobalVars,my_save_size,id_to_loc);
		Asynch_Set_Temp_Files(asynch,t_b,(void*)&t_b,0);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);

	if(my_rank == 0)	printf("\nTime for calculations: %f. Writing results to disk...\nCalculations complete! All done!\n",difftime(stop,start));

	//Clean up
	free(x_start);
	free(bool_above_gauges);
	free(above_gauges);
	free(d);
	free(d_full);
	for(i=0;i<numdata;i++)
		free(H[i]);
	free(H);
/*
	MatDestroy(HM);
	MatDestroy(HTH);
	VecDestroy(RHS);
	free(HM_els);
	free(cols_allstates);
	KSPDestroy(&ksp);
*/
	//for(i=0;i<numdata*max_steps_to_use;i++)
	//	free(HM_full[i]);
	free(HM_full);
	//for(i=0;i<numdata;i++)
	//	free(HM[i]);
	//free(HM);


	free(ids);
	free(data_locs);
	for(i=0;i<numlinks_obs;i++)
	{
		for(j=0;j<numsteps[i];j++)
			free(truesolution[i][j]);
		free(truesolution[i]);
	}
	free(truesolution);
	free(numsteps);

	//TAO Clean up
	VecRestoreArray(Guess,&minimizer);
	VecDestroy(&Lowerbds);
	VecDestroy(&Upperbds);
	VecDestroy(&Guess);
	TaoDestroy(&tao);
	PetscFinalize();

	//Asynch clean up
	Asynch_Delete_Temporary_Files(asynch);
	Asynch_Free(asynch);

	return 0;
}



//Evaluates sum (q_g(x_0) - d_full)^2, and its gradient.
PetscErrorCode EvaluateFunctionAndGradient(Tao tao, Vec X, PetscReal *f, Vec G, void *ptr)
{
	unsigned int i,j,k,m,n,l,counter;
	Link* current;
	double factor;
	time_t start,stop;
	time(&start);

	//Unpack ptr
	AppCtxFull* user = (AppCtxFull*) ptr;
	asynchsolver* asynch = user->asynch;
	unsigned int N = asynch->N;
	Link** sys = asynch->sys;
	UnivVars* GlobalVars = asynch->GlobalVars;
	double *HM_full = user->HM_full,*d = user->d,*x_start = user->x_start;		//!!!! HM_full here only needs to be 1 x allstates_needed !!!!
	unsigned int *above_gauges = user->above_gauges,num_above = user->num_above,*data_locs = user->data_locs;
	//unsigned int N = user->N,my_N = user->my_N,my_max_nodes = user->my_max_nodes;
	//UnivVars* GlobalVars = user->GlobalVars;
	//Link **sys = user->sys;
	//unsigned int *my_sys = user->my_sys,*assignments = user->assignments,**id_to_loc = user->id_to_loc;
	//TransData* my_data = user->my_data;
	unsigned int problem_dim = user->problem_dim,allstates = user->allstates,inc = user->inc,max_steps_to_use = user->max_steps_to_use,steps_to_use = user->steps_to_use,numdata = user->numdata;
	double t_b = user->t_b;
	//TempStorage* workspace = user->workspace;
	unsigned int allstates_needed = user->allstates_needed;

	//Init system
	*f = 0.0;
	double *x,*g,sum;
	VecGetArray(X,&x);
	VecGetArray(G,&g);
	for(k=0;k<allstates_needed;k++)	g[k] = 0.0;
	for(i=0;i<num_above;i++)
	{
		for(j=0;j<problem_dim;j++)
			x_start[above_gauges[i]*problem_dim+j] = x[i*problem_dim+j];
	}

if(print_out)
{
printf("Start %u\n",eval_counter++);
printf("x\n");
Print_VECTOR(x,allstates_needed);
printf("x_start\n");
Print_VECTOR(x_start,allstates);
}

	//Asynch_Set_System_State(asynch,t_b,VEC** backup);	//!!!! Should have something like this... !!!!
	ResetSysLS(sys,N,GlobalVars,t_b,x_start,problem_dim,GlobalVars->num_forcings);
	for(i=0;i<N;i++)				//!!!! Put this into ResetSysLS? !!!!
		ReadInitData(GlobalVars->global_params,sys[i]->params,sys[i]->iparams,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type,GlobalVars->diff_start,GlobalVars->no_ini_start,NULL);	//!!!! Very inefficient. Too many checks. !!!!
	for(i=0;i<asynch->GlobalVars->num_forcings;i++)
	{
		if(asynch->forcings[i]->flag == 3)
		{
			Asynch_Set_Forcing_State(asynch,i,t_b,asynch->forcings[i]->first_file,asynch->forcings[i]->last_file);
			//Asynch_Set_First_Rainfall_Timestamp(asynch,asynch->forcings[i]->first_file,i);
			//Asynch_Set_Last_Rainfall_Timestamp(asynch,asynch->forcings[i]->last_file,i);
		}
	}
/*
if(print_out)
{
printf("Start %u\n",eval_counter++);
Print_VECTOR(x,allstates);
for(i=0;i<N;i++)
	Print_Vector(sys[i]->list->tail->y_approx);
printf("*******\n");
}
*/

printf("HM_full\n");

	unsigned int max_or_steps = (steps_to_use == max_steps_to_use) ? max_steps_to_use : steps_to_use;
	for(i=0;i<max_or_steps;i++)	//!!!! Start at i=1? For i = 0, I don't think we need to set anything... !!!!		//HM_full here holds the values of M that are needed
	{
		GlobalVars->maxtime = t_b + (i) * inc;
		if(i)	Asynch_Advance(asynch,0);

		//Build M's and f
		for(j=0;j<numdata;j++)
		{
			current = sys[data_locs[j]];	//!!!! Assumes only discharges !!!!
			for(n=0;n<num_above;n++)
			{
				for(l=0;l<problem_dim;l++)
					HM_full[above_gauges[n]*problem_dim+l] = 0.0;
			}

			//f
			*f += sq(current->list->tail->y_approx->ve[0] - d[i*numdata + j]);

			//M
			//From my link
			for(n=0;n<problem_dim;n++)
				HM_full[current->location*problem_dim + n] = current->list->tail->y_approx->ve[problem_dim + (problem_dim-1) + n];
				//HM_full[i*numdata + j][current->location*problem_dim + n] = current->list->tail->y_approx->ve[problem_dim + (problem_dim-1) + n];

			//From parents
			counter = 3*problem_dim - 1;
			for(n=0;n<current->numparents;n++)
			{
				for(l=0;l<current->numupstream[n];l++)
				{
					for(m=0;m<problem_dim;m++)
						HM_full[current->upstream[n][l]*problem_dim + m] = current->list->tail->y_approx->ve[counter++];
						//HM_full[i*numdata + j][current->upstream[n][l]*problem_dim + m] = current->list->tail->y_approx->ve[counter++];
				}
			}
/*
if(print_out)
{
printf("q_%i(%e) = %.16e\n",current->ID,current->last_t,current->list->tail->y_approx->ve[0]);
}
*/

Print_VECTOR(HM_full,allstates_needed);

			//Build gradient
			factor = 2.0 * (current->list->tail->y_approx->ve[0] - d[i*numdata+j]);
			for(k=0;k<num_above;k++)
			{
				for(n=0;n<problem_dim;n++)
					g[k*problem_dim + n] += factor * HM_full[above_gauges[k]*problem_dim + n];
					//g[k*problem_dim + n] += 2.0 * (current->list->tail->y_approx->ve[0] - d[i*numdata+j]) * HM_full[i*numdata + j][above_gauges[k]*problem_dim + n];
			}
		}
	}
if(print_out)
{
printf("f = %e\n",*f);
//printf("HM\n");
//Print_MATRIX(HM_full,numdata*max_steps_to_use,allstates);
printf("gradient\n");
Print_VECTOR(g,allstates_needed);
printf("\n");
//getchar();
}
	VecRestoreArray(X,&x);
	VecRestoreArray(G,&g);

	time(&stop);
	printf("Total time for function eval: %.0f\n",difftime(stop,start));

	return 0;
}

/*
//Evaluates the hessian of sum (HM_full * x - d_full)^2
PetscErrorCode EvaluateHessian(Tao tao,Vec x,Mat *H,Mat *Hpre,MatStructure *flag,void *ptr)
{
	AppCtx* user = (AppCtx*) ptr;

	unsigned int i,k,l;
	double** A = (double**) user->A;
	unsigned int allstates = user->allstates;
	unsigned int steps_to_use = user->steps_to_use;
	double** hess = user->hess;
	PetscInt* idmn = user->idmn;

	for(k=0;k<allstates;k++)
	{
		for(l=k;l<allstates;l++)
		{
			hess[k][l] = 0.0;
			for(i=0;i<steps_to_use;i++)
				hess[k][l] += 2.0 * A[i][k] * A[i][l];
			hess[l][k] = hess[k][l];
		}
	}

	MatSetValues(*H,allstates,idmn,allstates,idmn,user->hess_array,INSERT_VALUES);
	MatAssemblyBegin(*H,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*H,MAT_FINAL_ASSEMBLY);

	return 0;
}
*/

void LoadBounds(Vec* Lowerbds,Vec* Upperbds,unsigned int size)
{
	unsigned int i;
	PetscReal *lowerbds,*upperbds;
	VecGetArray(*Lowerbds,&lowerbds);
	VecGetArray(*Upperbds,&upperbds);

	for(i=0;i<size;i++)
	{
		if(i%2 == 0) lowerbds[i] = 1e0;
		else	lowerbds[i] = 0.0;
		//lowerbds[i] = 0.0;
		upperbds[i] = PETSC_INFINITY;
	}

lowerbds[0] = 1e-4;
lowerbds[2] = 1e-4;

	VecRestoreArray(*Lowerbds,&lowerbds);
	VecRestoreArray(*Upperbds,&upperbds);
}


void LoadBounds_percentage(Vec* Lowerbds,Vec* Upperbds,double* state,unsigned int size,double percent)
{
	unsigned int i;
	PetscReal *lowerbds,*upperbds;
	VecGetArray(*Lowerbds,&lowerbds);
	VecGetArray(*Upperbds,&upperbds);

	for(i=0;i<size;i++)
	{
//		lowerbds[i] = state[i] * percent;
		if(i%2 == 0)	lowerbds[i] = state[i] * percent;
		else	lowerbds[i] = 0.0;
		upperbds[i] = PETSC_INFINITY;
	}

//printf("lower\n");
//Print_VECTOR(lowerbds,size);


	VecRestoreArray(*Lowerbds,&lowerbds);
	VecRestoreArray(*Upperbds,&upperbds);
}
/*
!!!! Get the bounds scaling with area. What's up with MPI_Finalize?? !!!!
//Creates bounds that are within percent of the gauge readings, scaled by area
void CreateBounds_Scale(Tao tao,asynchsolver* asynch,Vec* Lowerbds,Vec* Upperbds,double* gauges,unsigned int num_gauges,double percent)
{


	TaoSetVariableBounds(tao,Lowerbds,Upperbds);
}
*/

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
/*
	for(i=0;i<numdata;i++)
	{
		for(current = sys[gauges[i]]->c;current && !(*bool_above_gauges)[current->location];current = current->c)
		{
			(*bool_above_gauges)[current->location] = 1;
			num_above++;
		}
	}
*/

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
	Link** sys = asynch->sys;
	unsigned int area_idx = GlobalVars->area_idx,problem_dim = GlobalVars->problem_dim;

	unsigned int i,j,loc,*upstream = (unsigned int*) malloc(N*sizeof(unsigned int)),num_upstream,prev_loc,curr_loc;
	short int *locs_set = (short int*) calloc(N,sizeof(short int));
	double *locs_newq = (double*) malloc(N*sizeof(double));
	unsigned int *counter = (unsigned int*) malloc(N*sizeof(unsigned int));
	Link* current;
	double ratio;

	for(i=0;i<N;i++)
	{
		locs_newq[i] = -1.0;
	}

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
					locs_newq[curr_loc] = current->params->ve[area_idx] * locs_newq[loc] / sys[loc]->params->ve[area_idx];
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
				ratio = locs_newq[current->location] / current->params->ve[area_idx];
				for(j=0;j<num_upstream;j++)
				{
					locs_set[upstream[j]] = 1;
					locs_newq[upstream[j]] = sys[upstream[j]]->params->ve[area_idx] * ratio;
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
	return 0;
}




//Read into memory the times and discharges stored in a .dat file.
double*** DownloadGaugeReadings(unsigned int start_time,unsigned int stop_time,unsigned int** id_to_loc,unsigned int N,unsigned int* numlinks,unsigned int** ids,unsigned int** locs,unsigned int** numsteps)
{
	unsigned int i,j,k;
	double*** data;
	PGresult *res;
	char query[1028];

	//Hard coding!
	unsigned int outlet = 307864;
	ConnData* conninfo = CreateConnData("dbname=model_ifc host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
	ConnectPGDB(conninfo);

	//Get link ids of gauges
	sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
		SELECT DISTINCT A.link_id FROM env_pois_adv AS A, subbasin AS B, sensor_data AS C \
		WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND B.link_id = A.link_id AND B.link_id != 311903;",outlet);
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
		sprintf(query,"SELECT (extract('epoch' FROM date_trunc('minute', C.time2)) - %u)/60 AS unix_time, getdischarges((dist_bottom - measured_dist)*0.0328084::real,A.id)*0.0283168 AS q \
			FROM env_pois_adv AS A, sensor_data AS C \
			WHERE A.type = 4 AND C.bridge_id = (A.foreign_id1)::integer AND A.link_id = %u AND to_timestamp(%u) <= C.time2 AND C.time2 <= to_timestamp(%u) \
			ORDER BY unix_time;",start_time,(*ids)[i],start_time,stop_time);
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
/*
printf("ids (%u)\n",*numlinks);
for(i=0;i<*numlinks;i++)
	printf("%u ",(*ids)[i]);
printf("\nsteps (%u)\n",*numlinks);
for(i=0;i<*numlinks;i++)
	printf("%u ",(*numsteps)[i]);
printf("\n\n");

for(i=0;i<(*numsteps)[0];i++)
	printf("%f %f\n",data[0][i][0],data[0][i][1]);

printf("Download finished...\n");
getchar();
*/
	return data;
}

