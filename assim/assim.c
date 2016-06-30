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
#include "parmathmethods.h"

int flaggy;
int flaggy2;

int my_rank;
int np;

/*
double ExactType0(double t,VEC* y_0,VEC* params,VEC* global_params,double stddev);
double ExactType0_yconn(double t,VEC* y_0,VEC* y_01,VEC* y_02,VEC* params,VEC* params1,VEC* params2,VEC* global_params,double stddev);
double ExactTestType0(unsigned int location,double t,double* y_0,double invtau,VEC* global_params,double stddev);
double ExactTestType1(unsigned int location,double t,double* y_0,double invtau,VEC* global_params,double stddev);
double ExactLeafType1(double t,double t_0,VEC* y_0,double rain_value,VEC* params,VEC* global_params,double stddev);
*/

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

	flaggy = 0;
	flaggy2 = 0;

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
	if(GlobalVars == NULL)
	{
		MPI_Finalize();
		return 1;
	}
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
	outputfile = PrepareTempFiles(sys,N,assignments,GlobalVars,save_list,save_size,my_save_size,filename,id_to_loc);
/*
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
	}
*/
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
	unsigned int ensemble_size = 10;		//Number of samples in the ensemble
	double rho = 10.0;
	unsigned int all_states = N*dim;
	VEC* true_init = v_get(all_states);
	for(i=0;i<all_states;i+=dim)
	{
		if(assignments[i/dim] == my_rank)
		{
			for(j=0;j<dim;j++)
				true_init->ve[i+j] = sys[i/dim]->list->head->y_approx->ve[j];
		}
	}
	//unsigned int data_locs[] = {1,3,4};	//Locations of links with data	!!!! Make sure these are locations and not IDs !!!!
	unsigned int data_locs[] = {0};
	unsigned int numdata = sizeof(data_locs)/sizeof(unsigned int);	//Number of links for which data is available
	double t_f = GlobalVars->maxtime;	//The actual final time. This gets overwritten in GlobalVars
	GlobalVars->maxtime = data_inc;
	double t_b;			//The background time

	//!!!! Not sure if multiplying N and numdata by dim is appropriate !!!!
	MAT* H = m_get(numdata,all_states);
	VEC* Rinv = v_get(numdata);
	MAT* Y_b = m_get(numdata,ensemble_size);
	VEC* y_b = v_get(numdata);		//Innovation average
	VEC* x_b = v_get(all_states);			//Background ensemble average
	MAT* X_b = m_get(all_states,ensemble_size);
	MAT* X_a = m_get(all_states,ensemble_size);
	MAT* C = m_get(ensemble_size,numdata);
	MAT* Pinv = m_get(ensemble_size,ensemble_size);
	MAT* P = m_get(ensemble_size,ensemble_size);
	MAT* V = m_get(ensemble_size,ensemble_size);		//Othogonal eigenvectors of Pinv
	VEC* D = v_get(ensemble_size);		//Eigenvalues of Pinv
	unsigned int* isuppz = (unsigned int*) malloc(2*ensemble_size*sizeof(ensemble_size));
	VEC* temp = v_get(ensemble_size);
	VEC* temp2 = v_get(ensemble_size);
	MAT* W = m_get(ensemble_size,ensemble_size);
	double sqrt_ensemble_size_m1 = sqrt(ensemble_size - 1.0);
	VEC* w = v_get(ensemble_size);
	VEC* yo = v_get(numdata);
	VEC* d = v_get(numdata);	//Innovations
	//double* sys_state = (double*) malloc(all_states*sizeof(double));
	VEC* sys_state = v_get(all_states);

	//Barrier and stop clock
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time = difftime(stop,start);
	if(my_rank == 0)	printf("Finished reading files. Total time: %f\n",total_time);

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,my_N);
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Blacs inits
	int context,nprow,npcol,my_row,my_col,row_blocking,col_blocking;
	unsigned int globalcol,globalrow;
	double randnum1,randnum2;
	//nprow = npcol = (int) rint(pow(np,.5));
	FindNPRowCol(&nprow,&npcol);
	Cblacs_get(0,0,&context);
	Cblacs_gridinit(&context,"R",nprow,npcol);
	my_row = my_rank / npcol;
	my_col = my_rank % npcol;
	row_blocking = 2;
	col_blocking = 2;

	//Initialize matrices
	for(i=0;i<numdata;i++)	H->me[i][data_locs[i]*dim] = 1.0;
	for(i=0;i<numdata;i++)	Rinv->ve[i] = 1.0/sq(data_stddev);

	//Set the initial conditions for the background ensemble
	for(j=0;j<ensemble_size;j++)
	{
		for(i=0;i<all_states;i+=dim)
		{
			if(assignments[i/dim] == my_rank)
			{
				for(l=0;l<dim;l++)
				{
					randnum1 = (rand()%100000)/100000.0 + .00000001;
					randnum2 = ((rand()%100000)/100000.0 + .00000001);
					X_b->me[i+l][j] = max( true_init->ve[l] + init_stddev->ve[l] * sqrt(-2.0*log( randnum1 ))*cos(TWO_PI * randnum2) , 0.0 );
					//X_b->me[i+l][j] = max( true_init->ve[l] + init_stddev->ve[l] * sqrt(-2.0*log( (rand()%100000)/100000.0 + .00000001 ))*cos(TWO_PI * ((rand()%100000)/100000.0 + .00000001)) , 0.0 );
					//X_b->me[i+l][j] = true_init->ve[i+l] + j/100.0;
				//if(X_b->me[i][j] < 0.0)		X_b->me[i][j] = 1e-6;
				}
			}
		}
	}

	//Send the initial conditions (X_b) to each process
	if(np > 1)
	{
		for(i=0;i<all_states;i+=dim)
			MPI_Bcast(X_b->me[i],ensemble_size*dim,MPI_DOUBLE,assignments[i/dim],MPI_COMM_WORLD);

		for(i=0;i<all_states;i+=dim)
			MPI_Bcast(&(true_init->ve[i]),dim,MPI_DOUBLE,assignments[i/dim],MPI_COMM_WORLD);
	}

//SetMat(X_b);

/*
printf("Here's X_b\n");
Print_Matrix(X_b);
printf("\n");
*/

	//Create X_blocal. !!!! Can the previous step be avoided? !!!!
	int* descX_b = malloc(9*sizeof(int));
	descX_b[0] = 1;
	descX_b[1] = context;
	descX_b[2] = all_states;
	descX_b[3] = ensemble_size;
	descX_b[4] = row_blocking;
	descX_b[5] = col_blocking;
	descX_b[6] = 0;
	descX_b[7] = 0;
	descX_b[8] = numroc_((int*)&all_states,&(descX_b[4]),&my_row,&(descX_b[6]),&nprow);
	MAT* X_blocal = m_get(numroc_((int*)&ensemble_size,&(descX_b[5]),&my_col,&(descX_b[7]),&npcol), numroc_((int*)&all_states,&(descX_b[4]),&my_row,&(descX_b[6]),&nprow));
	MAT* X_alocal = m_get(numroc_((int*)&ensemble_size,&(descX_b[5]),&my_col,&(descX_b[7]),&npcol), numroc_((int*)&all_states,&(descX_b[4]),&my_row,&(descX_b[6]),&nprow));

	for(i=0;i<all_states;i++)	//!!!! Can this be done better? !!!!
	{
		for(j=0;j<ensemble_size;j++)
		{
			if( (int)((i/row_blocking) % nprow) == my_row && (int)((j/col_blocking) % npcol) == my_col )
				X_blocal->me[(j/(npcol*col_blocking))*col_blocking+j%col_blocking][(i/(nprow*row_blocking))*row_blocking+i%row_blocking] = X_b->me[i][j];
		}
	}

	//Create Hlocal. !!!! Can this be avoided? !!!!
	int* descH = malloc(9*sizeof(int));
	descH[0] = 1;
	descH[1] = context;
	descH[2] = numdata;
	descH[3] = all_states;
	descH[4] = row_blocking;
	descH[5] = col_blocking;
	descH[6] = 0;
	descH[7] = 0;
	descH[8] = numroc_((int*)&numdata,&(descH[4]),&my_row,&(descH[6]),&nprow);
	MAT* Hlocal = m_get(numroc_((int*)&all_states,&(descH[5]),&my_col,&(descH[7]),&npcol),numroc_((int*)&numdata,&(descH[4]),&my_row,&(descH[6]),&nprow));

	for(i=0;i<numdata;i++)	//!!!! Can this be done better? !!!!
	{
		for(j=0;j<all_states;j++)
		{
			if( (int)((i/row_blocking) % nprow) == my_row && (int)((j/col_blocking) % npcol) == my_col )
				Hlocal->me[(j/(npcol*col_blocking))*col_blocking+j%col_blocking][(i/(nprow*row_blocking))*row_blocking+i%row_blocking] = H->me[i][j];
		}
	}


	//Read in the "true" solution
	unsigned int numlinks,*ids,*numsteps;
	double*** truesolution = ReadSolution("TempData/testobservations.dat",&numlinks,&ids,&numsteps);

	//Calculate number of times to loop
	unsigned int iterations = rint((t_f - GlobalVars->t_0) / data_inc);
	t_b = GlobalVars->t_0;

	//Make additional allocation for parallelization
	int* desca = malloc(9*sizeof(int));
	int* descz = malloc(9*sizeof(int));
	desca[0] = descz[0] = 1;
	desca[1] = descz[1] = context;
	desca[2] = descz[2] = ensemble_size;
	desca[3] = descz[3] = ensemble_size;
	desca[4] = descz[4] = row_blocking;
	desca[5] = descz[5] = col_blocking;
	desca[6] = descz[6] = 0;	//I think procs are numbered starting at 0
	desca[7] = descz[7] = 0;
	desca[8] = numroc_((int*)&ensemble_size,&(desca[4]),&my_row,&(desca[6]),&nprow);
	descz[8] = numroc_((int*)&ensemble_size,&(descz[4]),&my_row,&(descz[6]),&nprow);

	int lld_a = numroc_((int*)&ensemble_size,&(desca[4]),&my_row,&(desca[6]),&nprow);
	int loc_cols = numroc_((int*)&ensemble_size,&(desca[5]),&my_col,&(desca[7]),&npcol);
	MAT* Pinvlocal = m_get(loc_cols,lld_a);
	MAT* Plocal = m_get(loc_cols,lld_a);
	MAT* templocal = m_get(loc_cols,lld_a);
	MAT* Vlocal = m_get(loc_cols,lld_a);
	MAT* Wlocal = m_get(loc_cols,lld_a);
	//VEC* work = v_get((2 + 25*N + N*N)*2);
	VEC* work = v_get((2 + 25*ensemble_size + ensemble_size*ensemble_size)*2);
	//int NP0 = numroc_( ensemble_size, desca[4], 0, 0, nprow );
	//int MQ0 = numroc_( (int) max( ensemble_size, desca[4] ), desca[4], 0, 0, npcol );
	//VEC* work = v_get( 2 + 5*ensemble_size + (int) rint(max( 18*(ensemble_size+1), NP0 * MQ0 + 2 * desca[4] * desca[4] )) + (2 + (int) ceil(ensemble_size/np))*(ensemble_size+1) );
	//IVEC* iwork = iv_get((14*N + nprow*npcol + 5)*2);
	IVEC* iwork = iv_get((14*ensemble_size + nprow*npcol + 5)*2);
	VEC* temp_numdata = v_get(numdata);
	VEC* temp_all_states = v_get(all_states);

	//Create Y_blocal
	int* descY_b = malloc(9*sizeof(int));
	descY_b[0] = 1;
	descY_b[1] = context;
	descY_b[2] = numdata;
	descY_b[3] = ensemble_size;
	descY_b[4] = row_blocking;
	descY_b[5] = col_blocking;
	descY_b[6] = 0;
	descY_b[7] = 0;
	descY_b[8] = numroc_((int*)&numdata,&(descY_b[4]),&my_row,&(descY_b[6]),&nprow);
	MAT* Y_blocal = m_get(numroc_((int*)&ensemble_size,&(descY_b[5]),&my_col,&(descY_b[7]),&npcol),numroc_((int*)&numdata,&(descY_b[4]),&my_row,&(descY_b[6]),&nprow));

	//Create Clocal
	int* descC = malloc(9*sizeof(int));
	descC[0] = 1;
	descC[1] = context;
	descC[2] = ensemble_size;
	descC[3] = numdata;
	descC[4] = row_blocking;
	descC[5] = col_blocking;
	descC[6] = 0;
	descC[7] = 0;
	descC[8] = numroc_(&(descC[2]),&(descC[4]),&my_row,&(descC[6]),&nprow);
	MAT* Clocal = m_get(numroc_(&(descC[3]),&(descC[5]),&my_col,&(descC[7]),&npcol), numroc_(&(descC[2]),&(descC[4]),&my_row,&(descC[6]),&nprow));

	//Create yolocal and dlocal
	int* descyod = malloc(9*sizeof(int));
	descyod[0] = 1;
	descyod[1] = context;
	descyod[2] = numdata;
	descyod[3] = 1;
	descyod[4] = row_blocking;
	descyod[5] = col_blocking;
	descyod[6] = 0;
	descyod[7] = 0;
	descyod[8] = numroc_(&(descyod[2]),&(descyod[4]),&my_row,&(descyod[6]),&nprow);
	VEC* yolocal = v_get(numroc_(&(descyod[2]),&(descyod[4]),&my_row,&(descyod[6]),&nprow));
	VEC* dlocal = v_get(numroc_(&(descyod[2]),&(descyod[4]),&my_row,&(descyod[6]),&nprow));

	//Create tempw2local
	int* descw2 = malloc(9*sizeof(int));
	descw2[0] = 1;
	descw2[1] = context;
	descw2[2] = ensemble_size;
	descw2[3] = 1;
	descw2[4] = row_blocking;
	descw2[5] = col_blocking;
	descw2[6] = 0;
	descw2[7] = 0;
	descw2[8] = numroc_(&(descw2[2]),&(descw2[4]),&my_row,&(descw2[6]),&nprow);
	VEC* temp2local = v_get(numroc_(&(descw2[2]),&(descw2[4]),&my_row,&(descw2[6]),&nprow));
	VEC* wlocal = v_get(numroc_(&(descw2[2]),&(descw2[4]),&my_row,&(descw2[6]),&nprow));

	//Main loop
	for(k=0;k<iterations;k++)
	{
		//Step 1 ------------------------------

		//Apply H to the background ensemble
		//!!!! This can probably be done better. But it depends upon how H is defined. !!!!
		//mm_mlt(H,X_b,Y_b);
		parmm_mlt(Hlocal,descH,X_blocal,descX_b,Y_blocal,descY_b);


		//Calculate Y_b ensemble average
/*
		for(i=0;i<numdata;i++)
		{
			y_b->ve[i] = 0.0;
			for(j=0;j<ensemble_size;j++)
				y_b->ve[i] += Y_b->me[i][j];
			y_b->ve[i] *= 1.0/ensemble_size;
		}
*/
		parensemble_average(Y_blocal,descY_b,y_b,my_row,nprow,ensemble_size,temp_numdata,1);

		//Calculate the matrix Y_b
/*
		for(i=0;i<numdata;i++)
			for(j=0;j<ensemble_size;j++)
				Y_b->me[i][j] -= y_b->ve[i];
*/

/*
if(k == 32)
{
	printf("\ny_b:\n");
	Print_Vector(y_b);
}
*/


		//Step 2 ------------------------------

/*
		//Calculate background ensemble average
		for(i=0;i<all_states;i++)
		{
			x_b->ve[i] = 0.0;
			for(j=0;j<ensemble_size;j++)
				x_b->ve[i] += X_b->me[i][j];
			x_b->ve[i] *= 1.0/ensemble_size;
		}

		//Calculate the matrix X_b
		for(i=0;i<all_states;i++)
			for(j=0;j<ensemble_size;j++)
				X_b->me[i][j] -= x_b->ve[i];
*/

		//Calculate background ensemble average
		parensemble_average(X_blocal,descX_b,x_b,my_row,nprow,ensemble_size,temp_all_states,1);

		//Step 4 ------------------------------
/*
printf("Y_b: %f\n",Rinv->ve[0]);
parPrint_Matrix(Y_blocal,descY_b,my_row,my_col,nprow,npcol);
getchar();
*/
		//!!!! Can Rinv just be stored as a scalar? !!!!

		//Compute C
		//mTdiag_mlt(Y_b,Rinv,C);
		//parmTdiag_mlt(Y_blocal,descY_b,Rinv,Clocal,my_row,nprow);
		parTranspose(Y_blocal,descY_b,Rinv->ve[0],Clocal,descC);	//!!!! This assumes Rinv = const * I !!!!

		//Step 5 ------------------------------

		//Build Pinv
		//mm_mlt(C,Y_b,Pinv);
		parmm_mlt(Clocal,descC,Y_blocal,descY_b,Pinvlocal,desca);

		//for(i=0;i<ensemble_size;i++)
		//	Pinv->me[i][i] += (ensemble_size-1.0)/rho;
		parmaeye_add(Pinvlocal,desca,(ensemble_size-1.0)/rho,my_row,my_col,nprow,npcol);

//printf("Pinv: %u\n",work->dim);
//parPrint_Matrix(Pinvlocal,desca,my_row,my_col,nprow,npcol);

/*
		//Load up the local arrays
		for(i=0;i<ensemble_size;i++)	//!!!! Can this be done better? !!!!
		{
			for(j=0;j<ensemble_size;j++)
			{
				if( ((i/row_blocking) % nprow) == my_row && ((j/col_blocking) % npcol) == my_col )
					Pinvlocal->me[(j/(npcol*col_blocking))*col_blocking+j%col_blocking][(i/(nprow*row_blocking))*row_blocking+i%row_blocking] = Pinv->me[i][j];
			}
		}
*/

		//Compute eigendecomposition of Pinv
		//Diagonalize(Pinv,V,D,isuppz);
		Parallel_Diagonalize(Pinvlocal,Vlocal,D,desca,descz,work,iwork);

		//Compute P	!!!! This may not need to be done. Only used in Step 7 !!!!
		//VDinvVT_mlt(V,D,P,temp);
		parVDinvVT_mlt(Vlocal,desca,D,Plocal,templocal,my_col,npcol);

		//Step 6 ------------------------------

		//Build W
		//VsqrtDinvVT_mlt(V,D,W,temp);
		//for(i=0;i<ensemble_size;i++)
		//	W->me[i][i] *= sqrt_ensemble_size_m1;
		parVsqrtDinvVT_mlt(Vlocal,desca,D,Wlocal,templocal,my_col,npcol);

		unsigned int locrow = Vlocal->m;
		unsigned int loccol = Vlocal->n;

		for(i=0;i<locrow;i++)
		{
			for(j=0;j<loccol;j++)
			{
				globalcol = my_col*col_blocking + (i/col_blocking)*(col_blocking*npcol) + i%col_blocking;
				globalrow = my_row*row_blocking + (j/row_blocking)*(row_blocking*nprow) + j%row_blocking;
				if(globalcol == globalrow)
					Wlocal->me[i][j] *= sqrt_ensemble_size_m1;
			}
		}

		//Step 7 ------------------------------		

/*
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

		//Compute innovations
		v_sub(yo,y_b,d,0);
*/


		if(my_col == 0)
		{
			//Grab the observations
			for(i=0;i<yolocal->dim;i++)
			{
				globalrow = my_row*row_blocking + (i/row_blocking)*(row_blocking*nprow) + i%row_blocking;
				current = sys[data_locs[globalrow]];
				yolocal->ve[i] = FindDischarge(truesolution,current->ID,t_b,numlinks,ids,numsteps,data_stddev);
			}

			//Compute innovations
			//parvglobal_add(yolocal,descyod,-1.0,y_b,dlocal,my_col,npcol);
			parvglobal_add(yolocal,descyod,-1.0,y_b,dlocal,my_row,nprow);
		}
		MPI_Barrier(MPI_COMM_WORLD);
/*
if(k > 40)
{
	if(my_rank == 0)
	{
		printf("k = %u\n",k);
	}

	printf("y_b\n");
	Print_Vector(y_b);

	printf("Innovations:\n");
	parPrint_Vector(dlocal,descyod,my_row,my_col,nprow,npcol);
	getchar();
}
*/

		//Compute w
		//mv_mlt(C,d,temp);
		//mv_mlt(P,temp,w);

		parmv_mlt(Clocal,descC,dlocal,descyod,temp2local,descw2);
		parmv_mlt(Plocal,desca,temp2local,descw2,wlocal,descw2);

		//Add w to columns of W
		//for(i=0;i<ensemble_size;i++)
		//	for(j=0;j<ensemble_size;j++)
		//			W->me[i][j] += w->ve[i];

		parapply_vector(Wlocal,desca,wlocal,descw2,my_row,nprow,temp,temp2);

		//Step 8 ------------------------------

		//Build X_a
		//mm_mlt(X_b,W,X_a);
		//for(i=0;i<all_states;i++)
		//	for(j=0;j<ensemble_size;j++)
		//		X_a->me[i][j] += x_b->ve[i];

		parmm_mlt(X_blocal,descX_b,Wlocal,desca,X_alocal,descX_b);
		parapply_globalvector(X_alocal,descX_b,x_b,my_row,nprow);

		//Advances ----------------------------

		//Advance the analysis ensemble
		for(i=0;i<ensemble_size;i++)
		{
			//for(j=0;j<all_states;j++)	sys_state->ve[j] = X_a->me[j][i];	//!!!! Consider storing X_a and others as their transposes !!!!
			parglobalize_row(X_alocal,descX_b,my_row,my_col,nprow,npcol,i,sys_state->ve);
			ApplyState(sys,N,GlobalVars,sys_state->ve,t_b);
			AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,0,outputfile);
			Flush_TransData(my_data);
			GetState(sys,my_sys,N,my_N,assignments,GlobalVars,sys_state->ve);
			//for(j=0;j<all_states;j++)	X_b->me[j][i] = sys_state-ve[j];	//!!!! Same here !!!!
			parlocallize_row(X_blocal,descX_b,my_row,my_col,nprow,npcol,i,sys_state->ve);
		}
/*
if(my_rank == 0)
	printf("\nHere's X_a\n");
parPrint_Matrix(X_alocal,descX_b,my_row,my_col,nprow,npcol);
*/
		//Create the analysis
/*
		for(i=0;i<all_states;i++)
		{
			sys_state[i] = 0.0;
			for(j=0;j<ensemble_size;j++)
				sys_state[i] += X_a->me[i][j];
			sys_state[i] *= 1.0/ensemble_size;
		}
		ApplyState(sys,N,GlobalVars,sys_state,t_b,workspace);
*/

		parensemble_average(X_alocal,descX_b,sys_state,my_row,nprow,ensemble_size,temp_all_states,0);
/*
if(my_rank == 0)
{
	printf("\nHere's sys_state\n");
	Print_Vector(sys_state);
}
MPI_Barrier(MPI_COMM_WORLD);
sleep(1);
*/

/*
if(k == 32)
{
	printf("k = %u, here's sys_state\n",k);
	Print_Vector(sys_state);
	printf("\n");
}
*/
		ApplyState(sys,N,GlobalVars,sys_state->ve,t_b);
/*
if(k == 32)
{
	printf("Here's what's in sys:\n");
	for(l=0;l<N;l++)
	{
		printf("ID = %u\n",sys[l]->ID);
		Print_Vector(sys[l]->list->head->y_approx);
	}
	getchar();
}
*/

/*
MPI_Barrier(MPI_COMM_WORLD);
if(my_rank == 0) printf("Ok, done...\n");
sleep(500);
*/
		//Write the analysis to disk
		for(i=0;i<my_N;i++)		//Output the initial values to the file
		{
			current = sys[my_sys[i]];
			if(current->print_time > 0.0 && current->next_save - t_b < 1e-13 && t_b - current->next_save < 1e-13)
			{
				(current->disk_iterations)++;

				fsetpos(outputfile,&(current->pos));
				fwrite(&(current->next_save),sizeof(double),1,outputfile);
				for(j=0;j<GlobalVars->num_print;j++)
					fwrite(&(current->list->head->y_approx->ve[GlobalVars->print_indices[j]]),sizeof(double),1,outputfile);
				fgetpos(outputfile,&(current->pos));

/*
				fprintf(outputfile,"%i %.12e ",current->ID,current->next_save);
				for(j=0;j<GlobalVars->num_print;j++)
					fprintf(outputfile,"%.12e ",current->list->head->y_approx->ve[GlobalVars->dense_indices[j]]);
				fprintf(outputfile,"\n");
*/

				current->next_save += current->print_time;
			}
		}

		//Advance the system
		AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,1,outputfile);
		Flush_TransData(my_data);

		//Continue to next assimilation time
		GlobalVars->maxtime += data_inc;
		t_b += data_inc;
		MPI_Barrier(MPI_COMM_WORLD);
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
	v_free(sys_state);
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
	total_time += Process_Data(sys,GlobalVars,N,save_list,save_size,my_save_size,id_to_loc,assignments);

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

/*
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
*/

