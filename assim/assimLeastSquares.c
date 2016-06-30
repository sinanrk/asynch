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
#include "tao.h"

void MM_mult(double** A,double** B,double** C,unsigned int m,unsigned int inner,unsigned int n);
void VECTOR_Copy(double* u,double* v,unsigned int dim);
void LoadBounds(Vec* Lowerbds,Vec* Upperbds,unsigned int size);
PetscErrorCode EvaluateFunction(TaoSolver tao,Vec X,PetscReal* f,void *ptr);
PetscErrorCode EvaluateGradient(TaoSolver tao, Vec X, Vec G, void *ptr);
PetscErrorCode EvaluateHessian(TaoSolver tao,Vec x,Mat *H,Mat *Hpre,MatStructure *flag,void *ptr);
void Print_MATRIX(double** A,unsigned int m,unsigned int n);
void Print_VECTOR(double* v,unsigned int dim);


typedef struct
{
	double** A;
	double* d;
	unsigned int allstates;
	unsigned int steps_to_use;
	double** hess;
	double* hess_array;
	PetscInt* idmn;
} AppCtx;

int my_rank;
int np;

int main(int argc,char **argv)
{
	//MPI Stuff
	char help[] = "Here is a help message.\n";
	TaoInitialize(&argc,&argv,NULL,help);
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
		TaoFinalize();
		return 1;
	}

	//Asynch solver init

	//Declare variables
	time_t start,stop;
	unsigned int *my_sys,my_N,my_max_nodes,N,nummethods,my_save_size,save_size,i,j,k,l,m,n,*save_list,peaksave_size;
	int *assignments;
	char rkdfilename[256];
	short int *getting = NULL;
	TransData* my_data = NULL;
	RKMethod** AllMethods;
	Link** sys;
	Link* current;
	unsigned int** id_to_loc;
	ErrorData* GlobalErrors;
	UnivVars* GlobalVars;
	TempStorage* workspace;
	//PGconn* conn = NULL;
	ConnData* conninfo = CreateConnData("dbname=rm_model host=s-iihr58.iihr.uiowa.edu user=rainfall_feed password=r!Ain2012");
	char additional[16];	//For output filename

	//Read in .gbl file
	GlobalVars = Read_Global_Data(argv[1],&GlobalErrors,conninfo,rkdfilename);
	if(GlobalVars == NULL)	return 1;

	//Read in remaining data from files
	sys = Create_River_System_parallel(rkdfilename,&N,&my_sys,&my_N,&my_max_nodes,&my_data,&assignments,&getting,&AllMethods,&nummethods,GlobalVars,GlobalErrors,&save_list,&my_save_size,&save_size,&peaksave_size,&id_to_loc,conninfo,&workspace);
	if(sys == NULL)		return 1;

	//Put together the output filename string
	char outputfilename[GlobalVars->string_size],filename[GlobalVars->string_size];
	sprintf(outputfilename,"%s%s",GlobalVars->results_folder,GlobalVars->identifier);
	if(GlobalVars->print_par_flag == 1)
	{
		for(i=0;i<GlobalVars->global_params->dim;i++)
		{
			sprintf(filename,"_%.4e",GlobalVars->global_params->ve[i]);
			strcat(outputfilename,filename);
		}
	}

/*
	#ifdef PRINTPEAKFLOW
		//Setup the output peakflows files
		FILE* peakfile = NULL;
		char peakfilename[GlobalVars->string_size];
		sprintf(peakfilename,"%s.pea",outputfilename);

		if(my_rank == 0)
		{
			peakfile = fopen(peakfilename,"w");
			if(peakfile == NULL)
			{
				printf("Error creating peak file %s\n",peakfilename);
				return 1;
			}

			fprintf(peakfile,"%i\n%i\n\n",peaksave_size,GlobalVars->type);
			fclose(peakfile);
		}
	#endif
*/

	//Set print_time to t_0
	for(i=0;i<my_N;i++)
	{
		current = sys[my_sys[i]];
		if(current->save_flag)
		{
			current->next_save = current->last_t;
			current->disk_iterations = 0;
		}
	}

	//Setup temporary output data file
	FILE* outputfile = PrepareTempFiles(sys,N,assignments,GlobalVars,save_list,save_size,my_save_size,filename,NULL,id_to_loc);

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


	//Least Squares init

	//Observations
	unsigned int numlinks_obs,*ids,*numsteps;
	unsigned int data_locs[] = {2};	//Locations of links with data	!!!! Make sure these are locations and not IDs !!!!
	double*** truesolution = ReadSolution("TempData/yobservations.dat",&numlinks_obs,&ids,&numsteps);

	//Initialize choices
	double t_b = 0.0;
	double t_f = 300.0;
	double inc = 5.0;
	unsigned int numdata = 1;
	unsigned int dim = 2;
	unsigned int allstates = 6;
	double x_b[6] = {10, 0, 10, 0, 10, 0};
	double analysis[allstates];
	//unsigned int data = 2*2+1;
	double step = .1;
	double** H = (double**) malloc(numdata*sizeof(double*));
	for(i=0;i<numdata;i++)
		H[i] = (double*) calloc(allstates,sizeof(double));
	H[0][4] = 1.0;
	unsigned int max_steps_to_use = 10;

	//Other initializations
	unsigned int numstepstotake,steps_to_use,problem_dim = GlobalVars->problem_dim;
	unsigned int iterations = round((t_f - t_b) / inc);
	double* d = calloc(numdata,sizeof(double));
	unsigned int total_steps = round((inc - t_b) / step);
	unsigned int length_t = round((t_f - t_b) / step) + 1;

	unsigned int start_idx = 0;
	double* d_full = calloc(numdata*max_steps_to_use,sizeof(double));
	double** HM_full = (double**) malloc(numdata*max_steps_to_use*sizeof(double*));
	for(i=0;i<numdata*max_steps_to_use;i++)
		HM_full[i] = (double*) calloc(allstates,sizeof(double));
	double** M = (double**) malloc(allstates*sizeof(double*));
	for(i=0;i<allstates;i++)
		M[i] = (double*) calloc(allstates,sizeof(double));
	double** HM = (double**) malloc(numdata*sizeof(double*));
	for(i=0;i<numdata;i++)
		HM[i] = (double*) malloc(allstates*sizeof(double));

	//Initialize HM_full
	for(i=0;i<numdata;i++)
		for(j=0;j<allstates;j++)	HM_full[i][j] = H[i][j];	// !!!! This could probably be better... !!!!

	//Problem specific declarations
	double* ans_0 = (double*) calloc((allstates+dim+dim-1),sizeof(double));
	double* ans_1 = (double*) calloc((allstates+dim+dim-1),sizeof(double));
	double* ans_2 = (double*) calloc((allstates+dim+dim-1),sizeof(double));
	double* use_state_0 = (double*) malloc(allstates*sizeof(double));
	double* use_state_1 = (double*) malloc(allstates*sizeof(double));
	double* use_state_2 = (double*) malloc(allstates*sizeof(double));

	//Prep TAO
	if(my_rank == 0)	printf("\nPrepping TAO...\n");

	AppCtx user;
	user.A = HM_full;
	user.d = d_full;
	user.allstates = allstates;
	user.hess_array = (double*) malloc(allstates*allstates*sizeof(double));
	user.hess = (double**) malloc(allstates*sizeof(double*));
	for(i=0;i<allstates;i++)	user.hess[i] = &(user.hess_array[i*allstates]);
	user.idmn = (PetscInt*) malloc(allstates*sizeof(PetscInt));
	for(i=0;i<allstates;i++)	user.idmn[i] = i;
	Vec Lowerbds,Upperbds,Guess;
	Mat Hess;
	PetscInt indices[allstates];
	double* minimizer;
	for(i=0;i<allstates;i++)	indices[i] = i;
	TaoSolver tao;
	TaoCreate(PETSC_COMM_SELF,&tao);
	TaoSetType(tao,"tao_tron");
	//TaoSetType(tao,"tao_gpcg");
	//TaoSetType(tao,"tao_bqpip");
	VecCreateSeq(MPI_COMM_SELF,allstates,&Lowerbds);
	VecCreateSeq(MPI_COMM_SELF,allstates,&Upperbds);
	VecCreateSeq(MPI_COMM_SELF,allstates,&Guess);
	MatCreateSeqDense(MPI_COMM_SELF,allstates,allstates,PETSC_NULL,&Hess);
	LoadBounds(&Lowerbds,&Upperbds,allstates);
	TaoSetVariableBounds(tao,Lowerbds,Upperbds);
	TaoSetInitialVector(tao,Guess);
	TaoSetObjectiveRoutine(tao,EvaluateFunction,(void*)&user);
	TaoSetGradientRoutine(tao,EvaluateGradient,(void*)&user);
	TaoSetHessianRoutine(tao,Hess,Hess,EvaluateHessian,(void*)&user);

	VecGetArray(Guess,&minimizer);	//Grab solution

	//TaoSetTolerances(tao,1e-14,1e-12,1e-10,1e-10,0.0);

	double* zeros = malloc(allstates*sizeof(double));
	for(i=0;i<allstates;i++)	zeros[i] = 0.0;

	//Compute data assim model
	if(my_rank == 0)	printf("\nMaking calculations...\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	for(k=0;k<iterations;k++)
	{
		printf("\n\n*************************************\n");
		printf("Iteration %u/%u. Background = %e.\n",k,iterations,t_b);

		//steps_to_use = (k+1 < max_steps_to_use+1) ? k+1 : max_steps_to_use;
		steps_to_use = (k < max_steps_to_use) ? k+1 : max_steps_to_use;

		//Calculate new observations
		for(i=0;i<numdata;i++)	//!!!! Assuming 1 observation for now !!!!
		{
			current = sys[data_locs[i]];
			//d[i] = FindDischarge(truesolution,current->ID,0.0 + (k+1)*inc,numlinks_obs,ids,numsteps,0.0);
			d[i] = FindDischarge(truesolution,current->ID,0.0 + (k)*inc,numlinks_obs,ids,numsteps,0.0);
		}

		if(k >= max_steps_to_use)
		{
			for(j=0;j<(max_steps_to_use-1)*(numdata);j++)
				d_full[j] = d_full[j+1];
			for(j=0;j<numdata;j++)
				d_full[(max_steps_to_use-1)*(numdata) + j] = d[j];
		}
		else
		{
			for(j=0;j<numdata;j++)
				d_full[(steps_to_use-1)*numdata + j] = d[j];
		}


//printf("Guess %u\n",steps_to_use);
//Print_VECTOR(x_b,allstates);
//printf("indices\n");
//for(i=0;i<allstates;i++)	printf("%i ",indices[i]);
//printf("\n");
printf("HM\n");
Print_MATRIX(HM_full,max_steps_to_use,allstates);
printf("\n");
printf("d_full\n");
Print_VECTOR(d_full,steps_to_use);
printf("\n");


		//Calculate analysis
		user.steps_to_use = steps_to_use;
		VecSetValues(Guess,allstates,indices,x_b,INSERT_VALUES);
		CHKERRQ(TaoSolve(tao));
		TaoView(tao,PETSC_VIEWER_STDOUT_SELF);
		for(i=0;i<allstates;i++)	analysis[i] = minimizer[i];	//!!!! Can't minimizer just be analysis? !!!!


//if(k == 10)
{
printf("x_b\n");
Print_VECTOR(x_b,allstates);
printf("\n");

printf("analysis t_b = %e k = %u start_idx = %u last assim time = %e\n",t_b,k,start_idx,t_b + 5.0*(steps_to_use-1));
Print_VECTOR(analysis,allstates);
//getchar();
}

		//Form the new transition matrix

		//Load analysis into sys, and reset variational equation state
		ResetSysLS(sys,N,GlobalVars,t_b,analysis,problem_dim);
		for(i=0;i<N;i++)				//!!!! Put this into ResetSysLS? !!!!
			ReadInitData(GlobalVars->global_params,sys[i]->params,sys[i]->iparams,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type);

		unsigned int max_or_steps_p1 = (steps_to_use == max_steps_to_use) ? max_steps_to_use : steps_to_use + 1;
		for(i=0;i<max_or_steps_p1;i++)	//!!!! Start at i=1? For i = 0, I don't think we need to set anything... !!!!
		{
			//Write init condition to disk
			for(j=0;j<my_N;j++)
			{
				current = sys[my_sys[j]];
				if(current->save_flag && fabs(current->last_t - current->next_save) < 1e-10 )
				{
//if(current->ID == 0)
//printf("Assim: %e %e\n",current->next_save,current->list->tail->y_approx->ve[1]);
					fsetpos(outputfile,&(current->pos));
					fwrite(&(current->next_save),sizeof(double),1,outputfile);
					for(l=0;l<GlobalVars->num_print;l++)
						fwrite(&(current->list->tail->y_approx->ve[GlobalVars->print_indices[l]]),sizeof(double),1,outputfile);
					fgetpos(outputfile,&(current->pos));
					current->next_save += current->print_time;
					(current->disk_iterations)++;
				}
			}

			//GlobalVars->maxtime = t_b + (i+1) * inc;
			GlobalVars->maxtime = t_b + (i) * inc;
			if(i)	AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,1,outputfile);

			//Extract M !!!! Try to build HM directly !!!!   !!!! Not in parallel !!!!   !!!! Do we really need to compute M? Or just at links with sensors? !!!!
			for(n=0;n<N;n++)
			{
				current = sys[n];

				//From my link
				for(j=0;j<(problem_dim-1);j++)
					M[n*problem_dim + 1 + j][n*problem_dim + 1 + j] = current->list->tail->y_approx->ve[problem_dim + j];

				for(j=0;j<problem_dim;j++)
					M[n*problem_dim][n*problem_dim + j] = current->list->tail->y_approx->ve[problem_dim + (problem_dim-1) + j];

				//From parents
				unsigned int counter = 3*problem_dim - 1;
				for(j=0;j<current->numparents;j++)
				{
					for(l=0;l<current->numupstream[j];l++)
					{
						for(m=0;m<problem_dim;m++)
						{
							M[n*problem_dim][current->upstream[j][l]*problem_dim + m] = current->list->tail->y_approx->ve[counter];
							counter++;
						}
					}
				}
			}

printf("**********\n");
printf("i = %u\n",i);
Print_MATRIX(M,allstates,allstates);
printf("**********\n");

			//Form HM
			MM_mult(H,M,HM,numdata,allstates,allstates);

			//Form HM_full
			for(j=0;j<numdata;j++)
			{
				for(l=0;l<allstates;l++)
					HM_full[i*numdata+j][l] = HM[j][l];
			}

			//Save next x_b			!!!! Not in parallel !!!!
			if(steps_to_use == max_steps_to_use)
			{
				if(i == 1)
				{
//printf("!!!! Uh, using time %e for the background time. !!!!\n",sys[0]->list->head->t);
					for(j=0;j<N;j++)
						for(l=0;l<dim;l++)
							x_b[j*dim+l] = sys[j]->list->tail->y_approx->ve[l];
							//x_b[j*dim+l] = sys[j]->list->head->y_approx->ve[l];
				}
			}
			else
			{
				if(i == 0)
				{
					VECTOR_Copy(analysis,x_b,allstates);
/*
					for(j=0;j<N;j++)
						for(l=0;l<dim;l++)
							x_b[j*dim+l] = sys[j]->list->tail->y_approx->ve[l];
							//x_b[j*dim+l] = sys[j]->list->head->y_approx->ve[l];
*/
				}
			}
		}

		//Make predictions
		GlobalVars->maxtime = t_f;
		AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,1,outputfile);

		//Go to next time
		if(k >= max_steps_to_use-1)
		{
			start_idx += total_steps;
			t_b += inc;
		}
		else
		{
			start_idx = 0;			// !!!! Or do nothing at all? !!!!
		}

		//Reset the temporary files
		printf("Creating output file...\n");
		sprintf(additional,"%u",k);	//Add the iteration number to the output files
		Process_Data(sys,GlobalVars,N,save_list,save_size,my_save_size,id_to_loc,assignments,NULL,additional);
		SetTempFiles(t_b,sys,N,outputfile,GlobalVars,my_save_size,id_to_loc);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);

	if(my_rank == 0)	printf("\nTime for calculations: %f. Writing results to disk...\n",difftime(stop,start));
	if(my_rank == 0)	printf("Calculations complete! All done!\n");

	//Clean up
	free(user.hess_array);
	free(user.hess);
	free(user.idmn);
	free(d);
	free(d_full);
	for(i=0;i<numdata;i++)
		free(H[i]);
	free(H);
	for(i=0;i<numdata*max_steps_to_use;i++)
		free(HM_full[i]);
	free(HM_full);
	for(i=0;i<allstates;i++)
		free(M[i]);
	free(M);
	for(i=0;i<numdata;i++)
		free(HM[i]);
	free(HM);
	free(zeros);

	free(ans_0);
	free(ans_1);
	free(ans_2);
	free(use_state_0);
	free(use_state_1);
	free(use_state_2);

	free(ids);
	for(i=0;i<numlinks_obs;i++)
	{
		for(j=0;j<numsteps[i];j++)
			free(truesolution[i][j]);
		free(truesolution[i]);
	}
	free(truesolution);
	free(numsteps);

	//Asynch clean up
	TransData_Free(my_data);
	ConnData_Free(conninfo);
	Destroy_ErrorData(GlobalErrors);
	Destroy_Workspace(workspace,GlobalVars->max_s,GlobalVars->max_parents);
	free(workspace);
	free(getting);
	if(GlobalVars->rain_flag == 4)	Destroy_RainData(sys[my_sys[0]]->rain);
	for(i=0;i<N;i++)	Destroy_Link(sys[i],GlobalVars->iter_limit,rkdfilename[0] != '\0',GlobalVars);
	free(sys);
	free(my_sys);
	free(assignments);
	for(i=0;i<nummethods;i++)	Destroy_RKMethod(AllMethods[i]);
	free(AllMethods);
	if(save_list != NULL)	free(save_list);
	for(i=0;i<N;i++)
		free(id_to_loc[i]);
	free(id_to_loc);
	Destroy_UnivVars(GlobalVars);
	for(i=0;i<N;i++)	v_free(backup[i]);
	free(backup);
	fclose(outputfile);

	//TAO Clean up
	VecRestoreArray(Guess,&minimizer);
	VecDestroy(&Lowerbds);
	VecDestroy(&Upperbds);
	VecDestroy(&Guess);
	TaoDestroy(&tao);
	MatDestroy(&Hess);
	TaoFinalize();

	return 0;
}

//Evaluates sum (HM_full * x - d_full)^2
//[analysis resnorm residual exitflag output] = lsqnonneg(HM_full,d_full,options);
PetscErrorCode EvaluateFunction(TaoSolver tao,Vec X,PetscReal* f,void *ptr)
{
	AppCtx* user = (AppCtx*) ptr;

	unsigned int i,j;
	double** A = (double**) user->A;
	double* d = (double*) user->d;
	unsigned int allstates = user->allstates;
	unsigned int steps_to_use = user->steps_to_use;
	double *x,sum;
	VecGetArray(X,&x);
	//VecGetArray(F,&f);
	*f = 0.0;

	for(i=0;i<steps_to_use;i++)
	{
		sum = -d[i];
		for(j=0;j<allstates;j++)
			sum += A[i][j] * x[j];
		*f += sum * sum;
	}

	VecRestoreArray(X,&x);

	return 0;
}

//Evaluates the gradient of sum (HM_full * x - d_full)^2
PetscErrorCode EvaluateGradient(TaoSolver tao, Vec X, Vec G, void *ptr)
{
	AppCtx* user = (AppCtx*) ptr;

	unsigned int i,j;
	double** A = (double**) user->A;
	double* d = (double*) user->d;
	unsigned int allstates = user->allstates;
	unsigned int steps_to_use = user->steps_to_use;

	double *x,*g,sum;
	VecGetArray(X,&x);
	VecGetArray(G,&g);

	for(i=0;i<allstates;i++)
		g[i] = 0.0;

	for(i=0;i<steps_to_use;i++)
	{
		sum = -d[i];
		for(j=0;j<allstates;j++)
			sum += A[i][j] * x[j];
		for(j=0;j<allstates;j++)
			g[j] += 2.0*A[i][j]*sum;
	}

	VecRestoreArray(X,&x);
	VecRestoreArray(G,&g);
	return 0;
}

//Evaluates the hessian of sum (HM_full * x - d_full)^2
PetscErrorCode EvaluateHessian(TaoSolver tao,Vec x,Mat *H,Mat *Hpre,MatStructure *flag,void *ptr)
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


void LoadBounds(Vec* Lowerbds,Vec* Upperbds,unsigned int size)
{
	unsigned int i;
	PetscReal *lowerbds,*upperbds;
	VecGetArray(*Lowerbds,&lowerbds);
	VecGetArray(*Upperbds,&upperbds);

	for(i=0;i<size;i++)
	{
		if(i%2 == 0) lowerbds[i] = 1e-12;
		else	lowerbds[i] = 0.0;
		//lowerbds[i] = 0.0;
		upperbds[i] = TAO_INFINITY;
	}

	VecRestoreArray(*Lowerbds,&lowerbds);
	VecRestoreArray(*Upperbds,&upperbds);
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

