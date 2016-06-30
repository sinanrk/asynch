#include "rkmethods.h"
#include "system.h"
#include "comm.h"
#include "riversys.h"
#include "processdata.h"
#include <time.h>
#include <libpq-fe.h>
#include "mpi.h"
#include "misc.h"
#include "rainfall.h"
#include "assimmethods.h"
#include "solvers.h"

int flaggy;
int flaggy2;

int my_rank;
int np;

double ExactType0(double t,VEC* y_0,VEC* params,VEC* global_params,double stddev);
double ExactType0_yconn(double t,VEC* y_0,VEC* y_01,VEC* y_02,VEC* params,VEC* params1,VEC* params2,VEC* global_params,double stddev);
double ExactTestType0(unsigned int location,double t,double* y_0,double invtau,VEC* global_params,double stddev);
double ExactTestType1(unsigned int location,double t,double* y_0,double invtau,VEC* global_params,double stddev);
double ExactLeafType1(double t,double t_0,VEC* y_0,double rain_value,VEC* params,VEC* global_params,double stddev);

int main(int argc,char* argv[])
{
	//Initialize MPI stuff
	MPI_Init(&argc,&argv);
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
		MPI_Finalize();
		return 1;
	}

	//Declare variables
	double total_time = 0.0;
	time_t start,stop;
	unsigned int *my_sys,my_N,my_max_nodes,N,nummethods,my_save_size,save_size,i,j,k,l,*save_list,peaksave_size;
	int *assignments;
	char filename[256],rkdfilename[256];
	FILE* outputfile = NULL;
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

	start = time(NULL);

	//Read in .gbl file
	GlobalVars = Read_Global_Data(argv[1],&GlobalErrors,conninfo,rkdfilename);
	if(GlobalVars == NULL)	return 1;
	GlobalVars->assim_flag = 1;

	//Read in remaining data from files
	sys = Create_River_System_parallel(rkdfilename,&N,&my_sys,&my_N,&my_max_nodes,&my_data,&assignments,&getting,&AllMethods,&nummethods,GlobalVars,GlobalErrors,&save_list,&my_save_size,&save_size,&peaksave_size,&id_to_loc,conninfo,&workspace);
	if(sys == NULL)		return 1;

	//Put together the output filename string
	char outputfilename[128];
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
		char peakfilename[256];
		sprintf(peakfilename,"%s.pea",outputfilename);

		if(my_rank == 0)
		{
			peakfile = fopen(peakfilename,"w");
			if(peakfile == NULL)
			{
				printf("Error creating peak file %s\n",peakfilename);
				return 1;
			}

			//fprintf(peakfile,"%i\n%i\n\n",N,GlobalVars->type);
			fclose(peakfile);
		}
	#endif
*/
	//Setup temporary output data file
	if(my_save_size > 0)
	{
		sprintf(filename,"%s%s%.3i.dat",GlobalVars->temp_folder,GlobalVars->identifier,my_rank);
		outputfile = fopen(filename,"w");
		if(outputfile == NULL)
		{
			printf("Error creating output file %s.\n",filename);
			return 1;
		}

		fprintf(outputfile,"%i\n%i\n\n",my_save_size,GlobalVars->type);
/*
		for(i=0;i<my_N;i++)		//Output the initial values to the file
		{
			if(sys[my_sys[i]]->print_time > 0.0)
			{
				fprintf(outputfile,"%i %.12e ",sys[my_sys[i]]->ID,sys[my_sys[i]]->last_t);
				//for(j=0;j<GlobalVars->dim;j++)
				for(j=0;j<GlobalVars->problem_dim;j++)
					fprintf(outputfile,"%.12e ",sys[my_sys[i]]->list->head->y_approx->ve[j]);
				fprintf(outputfile,"\n");
			}
		}
*/
	}

	//Set the seed
	srand(time(NULL));

	//Initial some assim specific variables
	double TWO_PI = 2.0 * 3.141592653589;
	unsigned int dim = GlobalVars->dim;
	double data_stddev = 1e-5;		//The standard deviation of the data
	VEC* init_stddev = v_get(dim);		//The standard deviation of the init conditions
	init_stddev->ve[0] = 1e-4;
	init_stddev->ve[1] = 1e-3;
	double data_inc = 5.0;			//Time increment when new data becomes available
	//unsigned int ensemble_size = 40;		//Number of samples in the ensemble
	unsigned int ensemble_size = 10;
	double rho = 10.0;
	VEC* true_init = v_get(N*dim);
	for(i=0;i<N*dim;i+=dim)
	{
		if(assignments[i/dim] == my_rank)
		{
			for(j=0;j<dim;j++)
				true_init->ve[i+j] = sys[i/dim]->list->head->y_approx->ve[j];
		}
	}
	VEC* data = workspace->temp;	//Memory to hold data
	unsigned int data_locs[] = {1,3,4};	//Locations of links with data	!!!! Make sure these are locations and not IDs !!!!
	//unsigned int data_locs[] = {1};
	unsigned int numdata = sizeof(data_locs)/sizeof(unsigned int);	//Number of links for which data is available
	double t_f = GlobalVars->maxtime;	//The actual final time. This gets overwritten in GlobalVars
	GlobalVars->maxtime = data_inc;
	double t_b;			//The background time
	double backup_time;

	//!!!! Not sure if multiplying N and numdata by dim is appropriate !!!!
	MAT* H = m_get(numdata,N*dim);
	VEC* Rinv = v_get(numdata);
	MAT* Y_b = m_get(numdata,ensemble_size);
	VEC* y_b = v_get(numdata);		//Innovation average
	VEC* x_b = v_get(N*dim);			//Background ensemble average
	MAT* X_b = m_get(N*dim,ensemble_size);
	MAT* X_a = m_get(N*dim,ensemble_size);
	MAT* C = m_get(ensemble_size,numdata);
	MAT* Pinv = m_get(ensemble_size,ensemble_size);
	MAT* P = m_get(ensemble_size,ensemble_size);
	MAT* V = m_get(ensemble_size,ensemble_size);		//Othogonal eigenvectors of Pinv
	VEC* D = v_get(ensemble_size);		//Eigenvalues of Pinv
	unsigned int* isuppz = (unsigned int*) malloc(2*ensemble_size*sizeof(ensemble_size));
	VEC* temp = v_get(ensemble_size);
	MAT* W = m_get(ensemble_size,ensemble_size);
	double sqrt_ensemble_size_m1 = sqrt(ensemble_size - 1.0);
	VEC* w = v_get(ensemble_size);
	VEC* yo = v_get(numdata);
	VEC* d = v_get(numdata);	//Innovations
	double* sys_state = (double*) malloc(N*dim*sizeof(double));

	//Barrier and stop clock
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time = difftime(stop,start);
	if(my_rank == 0)	printf("Finished reading files. Total time: %f\n",total_time);

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,my_N);
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Initialize matrices
	for(i=0;i<numdata;i++)	H->me[i][data_locs[i]*dim] = 1.0;
	for(i=0;i<numdata;i++)	Rinv->ve[i] = 1.0/sq(data_stddev);

	//Set the initial conditions for the background ensemble
	double randnum1,randnum2;
	for(j=0;j<ensemble_size;j++)
	{
		for(i=0;i<N*dim;i+=dim)
		{
			if(assignments[i/dim] == my_rank)
			{
				for(l=0;l<dim;l++)
				{
					randnum1 = (rand()%100000)/100000.0 + .00000001;
					randnum2 = ((rand()%100000)/100000.0 + .00000001);
					X_b->me[i+l][j] = max( true_init->ve[l] + init_stddev->ve[l] * sqrt(-2.0*log( randnum1 ))*cos(TWO_PI * randnum2) , 0.0 );
				}
					//X_b->me[i+l][j] = true_init->ve[i+l] + j/100.0;
				//if(X_b->me[i][j] < 0.0)		X_b->me[i][j] = 1e-6;
			}
		}
	}

	//Send the initial conditions (X_b) to each process
	if(np > 1)
	{
		for(i=0;i<N*dim;i+=dim)
			MPI_Bcast(X_b->me[i],ensemble_size*dim,MPI_DOUBLE,assignments[i/dim],MPI_COMM_WORLD);

		for(i=0;i<N*dim;i+=dim)
			MPI_Bcast(&(true_init->ve[i]),dim,MPI_DOUBLE,assignments[i/dim],MPI_COMM_WORLD);
	}

//Print_MatrixC(X_b);
SetMat(X_b);
//getchar();

	//Read in the "true" solution
	unsigned int numlinks,*ids,*numsteps;
	double*** truesolution = ReadSolution("TempData/testobservations.dat",&numlinks,&ids,&numsteps);

	//Calculate number of times to loop
	unsigned int iterations = rint((t_f - GlobalVars->t_0) / data_inc);
	t_b = GlobalVars->t_0;

int stophere = 0;

	//Main loop
	for(k=0;k<iterations;k++)
	{
		//Step 1 ------------------------------
/*
//if(k == 18)
{
printf("X_b   t_b = %f 7: %f  16: %f\n",t_b,X_b->me[2*5][7],X_b->me[2*5][16]);
}
*/
/*
if(k == stophere)
{
	printf("X_b %f\n",t_b);
	Print_Matrix(X_b);
	//getchar();
}
*/


		//Apply H to the background ensemble
		mm_mlt(H,X_b,Y_b);

		//Calculate Y_b ensemble average
		for(i=0;i<numdata;i++)
		{
			y_b->ve[i] = 0.0;
			for(j=0;j<ensemble_size;j++)
				y_b->ve[i] += Y_b->me[i][j];
			y_b->ve[i] *= 1.0/ensemble_size;
		}

		//Calculate the matrix Y_b
		for(i=0;i<numdata;i++)
			for(j=0;j<ensemble_size;j++)
				Y_b->me[i][j] -= y_b->ve[i];

		//Step 2 ------------------------------

		//Calculate background ensemble average
		for(i=0;i<N*dim;i++)
		{
			x_b->ve[i] = 0.0;
			for(j=0;j<ensemble_size;j++)
				x_b->ve[i] += X_b->me[i][j];
			x_b->ve[i] *= 1.0/ensemble_size;
		}

		//Calculate the matrix X_b
		for(i=0;i<N*dim;i++)
			for(j=0;j<ensemble_size;j++)
				X_b->me[i][j] -= x_b->ve[i];

		//Step 4 ------------------------------

		//Compute C
		mTdiag_mlt(Y_b,Rinv,C);

		//Step 5 ------------------------------

		//Build Pinv
		mm_mlt(C,Y_b,Pinv);

		for(i=0;i<ensemble_size;i++)
			Pinv->me[i][i] += (ensemble_size-1.0)/rho;
/*
if(k == stophere)
{
printf("Pinv is\n");
Print_Matrix(Pinv);
getchar();
}
*/
		//Compute eigendecomposition of Pinv
		Diagonalize(Pinv,V,D,isuppz);

		//Compute P	!!!! This may not need to be done. Only used in Step 7 !!!!
		VDinvVT_mlt(V,D,P,temp);

		//Step 6 ------------------------------

		//Build W
		VsqrtDinvVT_mlt(V,D,W,temp);
		for(i=0;i<ensemble_size;i++)
			W->me[i][i] *= sqrt_ensemble_size_m1;

		//Step 7 ------------------------------

		//Grab the observations
		for(i=0;i<numdata;i++)
		{
			current = sys[data_locs[i]];
			yo->ve[i] = FindDischarge(truesolution,current->ID,t_b,numlinks,ids,numsteps,data_stddev);
			//yo->ve[i] = ExactLeafType1(t_b,double t_0,VEC* y_0,double rain_value,current->params->params,GlobalVars->global_params,data_stddev);
			//yo->ve[i] = ExactTestType1(data_locs[i],t_b,true_init->ve,current->params->ve[12],GlobalVars->global_params,data_stddev);
			//yo->ve[i] = ExactTestType0(data_locs[i],t_b,true_init->ve,current->params->ve[12],GlobalVars->global_params,data_stddev);
			//yo->ve[i] = ExactType0(t_b,true_init,current->params,GlobalVars->global_params,data_stddev);
		}
/*
if(k == stophere)
{
printf("\nHere's yo %f\n",t_b);
Print_Vector(yo);
//getchar();
}
*/
		//Compute innovations
		v_sub(yo,y_b,d,0);

		//Compute w
		mv_mlt(C,d,temp);
		mv_mlt(P,temp,w);
/*
if(k == stophere)
{
printf("\nHere's W %f\n",t_b);
Print_Matrix(W);
printf("\nHere's w %f\n",t_b);
Print_Vector(w);
//getchar();
}
*/
		//Add w to columns of W
		for(i=0;i<ensemble_size;i++)
			for(j=0;j<ensemble_size;j++)
					W->me[i][j] += w->ve[i];

		//Step 8 ------------------------------
/*
if(k == 2)
{
printf("\n\nX_b is\n");
Print_Matrix(X_b);
printf("\n\nW is\n");
Print_Matrix(W);
}
*/
		//Build X_a
		mm_mlt(X_b,W,X_a);
/*
if(k == 2)
{
printf("\n\nX_a is\n");
Print_Matrix(X_a);
}
*/
		for(i=0;i<N*dim;i++)
			for(j=0;j<ensemble_size;j++)
				X_a->me[i][j] += x_b->ve[i];
/*
if(k == stophere)
{
printf("\nx_b is\n");
Print_Vector(x_b);
printf("\n\nX_a after is %f\n",t_b);
Print_Matrix(X_a);
getchar();
}
*/
		//Advances ----------------------------

		//Advance the analysis ensemble
		for(i=0;i<ensemble_size;i++)
		{
			for(j=0;j<N*dim;j++)	sys_state[j] = X_a->me[j][i];	//!!!! Consider storing X_a and others as their transposes !!!!
			ApplyState(sys,N,GlobalVars,sys_state,t_b,workspace);
			AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,0,outputfile);
			Flush_TransData2(my_data);
			GetState(sys,my_sys,N,my_N,assignments,GlobalVars,sys_state);
			for(j=0;j<N*dim;j++)	X_b->me[j][i] = sys_state[j];	//!!!! Same here !!!!
		}

		//Create the analysis
		for(i=0;i<N*dim;i++)
		{
			sys_state[i] = 0.0;
			for(j=0;j<ensemble_size;j++)
				sys_state[i] += X_a->me[i][j];
			sys_state[i] *= 1.0/ensemble_size;
		}
		ApplyState(sys,N,GlobalVars,sys_state,t_b,workspace);
/*
//if(k == 18)
{
printf("X_a   t_b = %f 7: %f  16: %f\n",t_b,X_a->me[2*5][7],X_a->me[2*5][16]);
getchar();
}
*/
/*
if(k == 18)
{
	printf("X_a %f\n",t_b);
	Print_Matrix(X_a);
	getchar();
}
*/
/*
if(k == 18)
{
	printf("sys_state %.16f ID = %u size = %u\n",sys_state[5*2],sys[5]->ID,ensemble_size);

	for(l=0;l<N*dim;l++)
		printf("%.16f\n",sys_state[l]);
	getchar();
}
*/
		//Write the analysis to disk
		for(i=0;i<my_N;i++)		//Output the initial values to the file
		{
			current = sys[my_sys[i]];
			if(current->print_time > 0.0 && current->next_save - t_b < 1e-13 && t_b - current->next_save < 1e-13)
			{
				(current->disk_iterations)++;
				fprintf(outputfile,"%i %.12e ",current->ID,current->next_save);
				//for(j=0;j<GlobalVars->dim;j++)
				//	fprintf(outputfile,"%.12e ",current->list->tail->y_approx->ve[j]);
				for(j=0;j<GlobalVars->num_print;j++)
					fprintf(outputfile,"%.12e ",current->list->head->y_approx->ve[GlobalVars->dense_indices[j]]);
				fprintf(outputfile,"\n");

				current->next_save += current->print_time;
			}
		}

		//Advance the system
		AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,1,outputfile);
		Flush_TransData2(my_data);

		//Continue to next assimilation time
		GlobalVars->maxtime += data_inc;
		t_b += data_inc;
		MPI_Barrier(MPI_COMM_WORLD);
/*
for(i=0;i<X_a->m;i++)
{
	for(j=0;j<X_a->n;j++)
		if(X_a->me[i][j] < 0.0)	X_a->me[i][j] = 0.0;
}
*/
	}

	//Free some memory
	MPI_Barrier(MPI_COMM_WORLD);
	TransData_Free(my_data);
	m_free(H);
	v_free(true_init);
	v_free(Rinv);
	m_free(Y_b);
	v_free(y_b);
	m_free(X_b);
	v_free(x_b);
	m_free(X_a);
	m_free(C);
	m_free(Pinv);
	m_free(P);
	m_free(V);
	v_free(D);
	free(isuppz);
	m_free(W);
	v_free(w);
	v_free(yo);
	v_free(d);
	free(sys_state);
	v_free(temp);
	v_free(init_stddev);
	for(i=0;i<numlinks;i++)
	{
		for(j=0;j<numsteps[i];j++)
			free(truesolution[i][j]);
		free(truesolution[i]);
	}
	free(truesolution);
	free(ids);
	free(numsteps);


/*
	#ifdef PRINTPEAKFLOW
		//Write number of iterations to the peak data to peakfile
		for(i=0;i<np;i++)
		{
			if(my_rank == i)
			{
				peakfile = fopen(peakfilename,"a");
				for(j=0;j<my_N;j++)
				{
					current = sys[my_sys[j]];
					if(current->peak_flag)
						fprintf(peakfile,"%i %.4f %.8f %.8f\n",current->ID,current->params->ve[0],current->peak_time,current->peak_value->ve[0]);
				}
				fclose(peakfile);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	#endif
*/
	//Stop the clock
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time += difftime(stop,start);

	//Output some data
	if(sys[my_sys[0]]->c == NULL)
	{
		printf("%i: The answer at ID %i at time %.12f is\n",my_rank,sys[my_sys[0]]->ID,sys[my_sys[0]]->last_t);
		Print_Vector(sys[my_sys[0]]->list->tail->y_approx);
		printf("Total time for calculations: %f\n",difftime(stop,start));
	}

	//Cleanup
	ConnData_Free(conninfo);
	Destroy_ErrorData(GlobalErrors);
	Destroy_Workspace(workspace,GlobalVars->max_s,GlobalVars->max_parents);
	free(workspace);
	free(getting);
	if(my_save_size > 0)	fclose(outputfile);
	MPI_Barrier(MPI_COMM_WORLD);	//Very important: make sure all processes have flushed their buffers before data is processed.

	//Process the data
	total_time += Process_Data(sys,GlobalVars,my_sys,N,my_N,save_list,save_size,my_save_size,id_to_loc,assignments);

	fflush(stdout);
	if(my_rank == 0)	printf("The total time for the entire program: %f\n\n",total_time);

	//Last bit of cleanup
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

	MPI_Finalize();
	return 0;
}

//For leaves
double ExactType0(double t,VEC* y_0,VEC* params,VEC* global_params,double stddev)
{
	double lambda_1 = global_params->ve[1];
	double invtau = params->ve[12];
	double TWO_PI = 2.0 * 3.141592653589;
	double mean;
	if(lambda_1 > 1e-12)
		mean = pow(1.0 + lambda_1 * invtau * pow(y_0->ve[0],lambda_1) * t,-1.0/lambda_1) * y_0->ve[0];
	else
		mean = y_0->ve[0] * exp(-invtau * t);

	return mean + stddev*sqrt(-2.0*log( (rand()%100000)/100000.0 + .00000001 ))*cos(TWO_PI * ((rand()%100000)/100000.0 + .00000001));
}

//Assumes lambda_1 = lambda_2 = 0.0
double ExactType0_yconn(double t,VEC* y_0,VEC* y_01,VEC* y_02,VEC* params,VEC* params1,VEC* params2,VEC* global_params,double stddev)
{
	double invtau = params->ve[12];
	double invtau1 = params1->ve[12];
	double invtau2 = params2->ve[12];
	double TWO_PI = 2.0 * 3.141592653589;
	double k1 = 1.0/(1.0 - invtau1 * 1.0/invtau);
	double k2 = 1.0/(1.0 - invtau2 * 1.0/invtau);
	double mean;

	if(fabs(invtau - invtau1) < 1e-12 && fabs(invtau1 - invtau2) < 1e-12)
		mean = invtau * (y_01->ve[0] + y_02->ve[0]) * t * exp(-invtau*t) + y_0->ve[0] * exp(-invtau*t);
	else
		mean = y_01->ve[0] * k1 * exp(-invtau1*t) + y_02->ve[0] * k2 * exp(-invtau2*t) + (y_0->ve[0] - k1*y_01->ve[0] - k2*y_02->ve[0])*exp(-invtau*t);

	return mean + stddev*sqrt(-2.0*log( (rand()%100000)/100000.0 + .00000001 ))*cos(TWO_PI * ((rand()%100000)/100000.0 + .00000001));
}

//This is designed to work with the test basin.
//Assumes lambda_1 = lambda_2 = 0.0 and all channel lengths are the same.
double ExactTestType0(unsigned int location,double t,double* y_0,double invtau,VEC* global_params,double stddev)
{
	unsigned int i;
	double TWO_PI = 2.0 * 3.141592653589;
	double mean;
	if(location == 1)
		mean = (1.0/6.0 * invtau*invtau*invtau * y_0[7] * t*t*t + .5 * invtau*invtau * (y_0[4] + y_0[6] + y_0[9] + y_0[10]) * t*t + invtau * (y_0[2] + y_0[3] + y_0[8]) * t + y_0[1]) * exp(-invtau*t);
	else if(location == 3)
		mean = (invtau * (y_0[9] + y_0[10]) * t + y_0[3]) * exp(-invtau*t);
	else if(location == 4)
		mean = (invtau * y_0[7]*t + y_0[4]) * exp(-invtau*t);
	else
	{
		printf("Error: invalid location. %u\n",location);
		return 0.0;
	}

	return mean + stddev*sqrt(-2.0*log( (rand()%100000)/100000.0 + .00000001 ))*cos(TWO_PI * ((rand()%100000)/100000.0 + .00000001));
}

//This is designed to work with the test basin.
//Assumes lambda_1 = lambda_2 = 0.0 and all channel lengths are the same.
double ExactTestType1(unsigned int location,double t,double* y_0,double invtau,VEC* global_params,double stddev)
{
	unsigned int i;
	double TWO_PI = 2.0 * 3.141592653589;
	double mean;
	if(location == 1)
		mean = (1.0/6.0 * invtau*invtau*invtau * y_0[7*2] * t*t*t + .5 * invtau*invtau * (y_0[4*2] + y_0[6*2] + y_0[9*2] + y_0[10*2]) * t*t + invtau * (y_0[2*2] + y_0[3*2] + y_0[8*2]) * t + y_0[1*2]) * exp(-invtau*t);
	else if(location == 3)
		mean = (invtau * (y_0[9*2] + y_0[10*2]) * t + y_0[3*2]) * exp(-invtau*t);
	else if(location == 4)
		mean = (invtau * y_0[7*2]*t + y_0[4*2]) * exp(-invtau*t);
	else
	{
		printf("Error: invalid location. %u\n",location);
		return 0.0;
	}

	return mean + stddev*sqrt(-2.0*log( (rand()%100000)/100000.0 + .00000001 ))*cos(TWO_PI * ((rand()%100000)/100000.0 + .00000001));
}

//Assumes lambda_1 = lambda_2 = 0.0 and all channel lengths are the same.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
double ExactLeafType1(double t,double t_0,VEC* y_0,double rain_value,VEC* params,VEC* global_params,double stddev)
{
	//!!!! Make this just read in values from disk !!!!

	double TWO_PI = 2.0 * 3.141592653589;

	double q_0 = y_0->ve[0];
	double s_0 = y_0->ve[1];

	double invtau = params->ve[12];
	double c_1 = params->ve[14];
	double c_3 = params->ve[16];
	double c_4 = params->ve[17];

	double expc_4 = exp(-c_4 *(t - t_0));
	double exptau = exp(-invtau *(t - t_0));
	double alpha = c_1 * (s_0 - c_3/c_4 * rain_value) * expc_4;
	double beta = c_1 * c_3/c_4 * rain_value;
	double frac = 1.0/(invtau - c_4);

	double mean = frac * invtau * alpha * expc_4 + beta + (q_0 - frac*invtau*alpha*exp(c_4*t_0) - beta) * exptau;

	return mean + stddev*sqrt(-2.0*log( (rand()%100000)/100000.0 + .00000001 ))*cos(TWO_PI * ((rand()%100000)/100000.0 + .00000001));
}


