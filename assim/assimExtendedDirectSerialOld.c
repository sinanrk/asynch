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
#include "cblas.h"
#include "clapack.h"

int my_rank;
int np;

//double ExactType0(double t,VEC* y_0,VEC* params,VEC* global_params,double stddev);
//double ExactType0_yconn(double t,VEC* y_0,VEC* y_01,VEC* y_02,VEC* params,VEC* params1,VEC* params2,VEC* global_params,double stddev);
//double ExactTestType0(unsigned int location,double t,double* y_0,double invtau,VEC* global_params,double stddev);

VEC* dump;

void ExactSolver(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,unsigned int my_max_nodes,UnivVars* GlobalVars,int* assignments,unsigned int** id_to_loc,TempStorage* workspace,ConnData* conninfo,TransData* my_data,short int print_flag,FILE* outputfile);
void ExactYConn(double t,double t_0,VEC* y_0_0,VEC* y_0_1,VEC* y_0_2,VEC* params_0,VEC* params_1,VEC* params_2,VEC* global_params,double precip,VEC* ans_0,VEC* ans_1,VEC* ans_2);
double FindDischargeExact(double t,Link** sys,UnivVars* GlobalVars);

int flaggy;

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
	unsigned int *my_sys,my_N,my_max_nodes,N,nummethods,my_save_size,save_size,i,j,k,l,m,*save_list,peaksave_size;
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
	outputfile = PrepareTempFiles(sys,N,assignments,GlobalVars,save_list,save_size,my_save_size,filename,NULL,id_to_loc);
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
	unsigned int problem_dim = GlobalVars->problem_dim;
	double data_stddev = 1e-2;		//The standard deviation of the data
	//double data_stddev = 1e-8;
	VEC* init_stddev = v_get(problem_dim);		//The standard deviation of the init conditions
	init_stddev->ve[0] = 1e+0;
	init_stddev->ve[1] = 1e+0;
	//init_stddev->ve[0] = 1e+4;
	//init_stddev->ve[1] = 1e+3;
	double data_inc = 5.0;		//Time increment when new data becomes available
	unsigned int data_locs[] = {2};	//Locations of links with data	!!!! Make sure these are locations and not IDs !!!!
	//unsigned int data_locs[] = {3};
	//unsigned int data_locs[] = {1,3,4};
	unsigned int numdata = sizeof(data_locs)/sizeof(unsigned int);	//Number of links for which data is available
	double t_f = GlobalVars->maxtime;	//The actual final time. This gets overwritten in GlobalVars
	GlobalVars->maxtime = data_inc;
	double t_b;			//The background time
	unsigned int all_states = N*problem_dim;
	MAT* H = m_get(numdata,all_states);
	MAT* B = m_get(all_states,all_states);
	MAT* R = m_get(numdata,numdata);
	VEC* d = v_get(numdata);	//Innovations
	VEC* analysis = v_get(all_states);	//Memory to hold the analysis
	MAT* K_temp = m_get(numdata,numdata);
	MAT* K = m_get(all_states,numdata);
	MAT* M = m_get(all_states,all_states);
	MAT* A = m_get(all_states,all_states);
	MAT* BHT = m_get(all_states,numdata);
	MAT* temp_dd1 = m_get(numdata,numdata);
	MAT* temp_dd2 = m_get(numdata,numdata);
	MAT* temp_NN = m_get(all_states,all_states);
	VEC* y_b = v_get(all_states);
	int* ipiv = (int*) malloc(numdata*sizeof(int));
	double* y_0 = (double*) malloc(all_states*sizeof(double));
	VEC* true_init = v_get(all_states);
	for(i=0;i<all_states;i+=problem_dim)
	{
		if(assignments[i/problem_dim] == my_rank)
		{
			for(j=0;j<problem_dim;j++)
				true_init->ve[i+j] = sys[i/problem_dim]->list->head->y_approx->ve[j];
		}
	}

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
	for(i=0;i<numdata;i++)	H->me[i][data_locs[i]*problem_dim] = 1.0;
	for(i=0;i<all_states;i+=problem_dim)
		for(l=0;l<problem_dim;l++)	B->me[i+l][i+l] = sq(init_stddev->ve[l]) * 1e3;
	for(i=0;i<numdata;i++)	R->me[i][i] = sq(data_stddev);
	for(i=0;i<all_states;i++)	M->me[i][i] = 1.0;

	//Set init conditions
	double TWO_PI = 2.0 * 3.141592653589;
	double randnum1,randnum2;
	//for(i=0;i<my_N;i++)
	//	sys[my_sys[i]]->list->head->y_approx->ve[0] = true_init->ve[0] + init_stddev * sqrt(-2.0*log( (rand()%100000)/100000.0 + .00000001 ))*cos(TWO_PI * ((rand()%100000)/100000.0 + .00000001));
	for(i=0;i<my_N*problem_dim;i+=problem_dim)
	{
		for(l=0;l<problem_dim;l++)
		{
			randnum1 = (rand()%100000)/100000.0 + .00000001;
			randnum2 = ((rand()%100000)/100000.0 + .00000001);
			sys[my_sys[i/problem_dim]]->list->head->y_approx->ve[l] = max( true_init->ve[problem_dim*sys[my_sys[i/problem_dim]]->location + l] + init_stddev->ve[l] * sqrt(-2.0*log( randnum1 ))*cos(TWO_PI * randnum2) , 0.0 );
		}
	}

//double array[] = { 9.6679414435823290, 0.1323403653890748, 6.6969472914860830, 0.0000000000000000, 3.6547434282435329, 0.0000000000000000 };
//double array[] = { 9.6679414435823290, 0.0, 6.6969472914860830, 0.0000000000000000, 3.6547434282435329, 0.0000000000000000 };
//double array[] = { 10.0, 0.0, 10.0, 0.0, 10.0, 0.0 };
double array[N*problem_dim];
for(i=0;i<N;i++)
{
	array[2*i] = 10.0;
	//array[2*i] = 9.0;
	array[2*i+1] = 0.0;
}
	
for(i=0;i<N;i++)
	for(j=0;j<problem_dim;j++)
		sys[i]->list->head->y_approx->ve[j] = array[i*problem_dim+j];
/*
//Set B
for(i=0;i<all_states;i+=2)
{
	for(j=0;j<all_states;j+=2)
		B->me[i][j] = 1.0;
}
*/

/*
printf("{ ");
for(i=0;i<N;i++)
	printf("%.16f, %.16f, ",sys[i]->list->head->y_approx->ve[0],sys[i]->list->head->y_approx->ve[1]);
printf(" }\n");
*/
/*
for(i=0;i<N;i++)
{
	printf("ID = %u\n",sys[i]->ID);
	Print_Vector(sys[i]->list->head->y_approx);
}
getchar();
*/
	//Read in the "true" solution
	unsigned int numlinks,*ids,*numsteps;
	//double*** truesolution = ReadSolution("TempData/testobservations.dat",&numlinks,&ids,&numsteps);
	double*** truesolution = ReadSolution("TempData/yobservations.dat",&numlinks,&ids,&numsteps);

	//Calculate number of times to loop
	unsigned int iterations = rint((t_f - GlobalVars->t_0) / data_inc);
	t_b = GlobalVars->t_0;

int stophere = -1;
flaggy = -1;

for(i=0;i<N;i++)
{
	printf("%i\n",sys[i]->ID);
	Print_Vector(sys[i]->params);
}
getchar();
//for(i=0;i<N;i++)
//printf("%u: alpha = %.16e beta = %.16e gamma = %.16e invtau = %.16e\n",sys[i]->ID,sys[i]->params->ve[14],sys[i]->params->ve[16],sys[i]->params->ve[17],sys[i]->params->ve[12]);


	//Main loop
	//for(t_b = GlobalVars->t_0; t_b <= t_f; t_b += data_inc)
	for(k=0;k<iterations;k++)
	{
		//Set background
		for(i=0;i<N;i++)
			for(j=0;j<problem_dim;j++)
				//y_b->ve[i*problem_dim+j] = sys[i]->list->head->y_approx->ve[j];
				y_b->ve[i*problem_dim+j] = sys[i]->list->tail->y_approx->ve[j];

		//Calculate innovations
		for(i=0;i<numdata;i++)
		{
			current = sys[data_locs[i]];
			//d->ve[i] = FindDischarge(truesolution,current->ID,t_b,numlinks,ids,numsteps,data_stddev) - current->list->tail->y_approx->ve[0];
			d->ve[i] = FindDischargeExact(t_b,sys,GlobalVars) - current->list->tail->y_approx->ve[0];

			//d->ve[i] = ExactType0(t_b,true_init,current->params,GlobalVars->global_params,data_stddev) - current->list->tail->y_approx->ve[0];
			//d->ve[0] = ExactType0_yconn(t_b,true_init,true_init,true_init,current->params,current->parents[0]->params,current->parents[1]->params,GlobalVars->global_params,data_stddev) - current->list->tail->y_approx->ve[0];
			//d->ve[i] = ExactTestType0(data_locs[i],t_b,y_0,sys[my_sys[0]]->params->ve[12],GlobalVars->global_params,data_stddev) - current->list->tail->y_approx->ve[0];
/*
if( (k == stophere || stophere == -1) )
{
printf("Approximation at link %u: %.16f\n",current->ID,current->list->tail->y_approx->ve[0]);
}
*/
		}

/*
//Make B symmetric
for(i=0;i<all_states;i++)
{
	for(j=i+1;j<all_states;j++)
		B->me[j][i] = B->me[i][j];
}
*/

//d->ve[0] = 0.0;

if( (k == stophere || stophere == -1) )
{
printf("\n\n************************************\n");
printf("t_b = %f k = %u\n",t_b,k);
printf("d is\n");
Print_Vector(d);
//printf("R is\n");
//Print_Matrix(R);
printf("B is\n");
Print_Matrix(B);
printf("y_b is\n");
Print_Vector(y_b);
//printf("H is\n");
//Print_Matrix(H);
//getchar();
}


		//Calculate gain matrix
		mmT_mlt(B,H,BHT);			// !!!! Simplify because of H? !!!!
		mm_mlt(H,BHT,temp_dd2);		// !!!! Simplify because of H? !!!!
		for(i=0;i<R->m;i++)	temp_dd2->me[i][i] += R->me[i][i];
if( (k == stophere || stophere == -1) )
{
printf("Inverting\n");
Print_Matrix(temp_dd2);
}
		clapack_dgetrf(CblasRowMajor,numdata,numdata,temp_dd2->array,temp_dd2->m,ipiv);
		clapack_dgetri(CblasRowMajor,numdata,temp_dd2->array,temp_dd2->m,ipiv);
if( (k == stophere || stophere == -1) )
{
printf("Got\n");
Print_Matrix(temp_dd2);
}
		mm_mlt(BHT,temp_dd2,K);		//!!!! Need to use ipiv? !!!!
/*
if(K->me[5][0] < 0.0)
{
	K->me[5][0] *= -1;
	printf("Fired!!!\n");
}
*/

if( (k == stophere || stophere == -1) )
{
printf("K is\n");
Print_Matrix(K);
//getchar();
}

		//Calculate analysis
		mv_mlt(K,d,analysis);
		v_add(y_b,analysis,analysis,0);		//!!!! Can this be rewritten to solve a linear system? !!!!

if( (k == stophere || stophere == -1) )
{
printf("Analysis\n");
Print_Vector(analysis);
//getchar();
}

/*
if( (k == stophere || stophere == -1) && t_b > 59.999)
{
printf("analysis\n");
Print_Vector(analysis);
//getchar();
}
*/

		//Set analysis
		for(i=0;i<N;i++)
		{
			for(j=0;j<problem_dim;j++)
				sys[i]->list->head->y_approx->ve[j] = analysis->ve[i*problem_dim+j];
//			CheckConsistency(sys[i]->list->head->y_approx,GlobalVars->dim,GlobalVars->type,sys[i]->params,GlobalVars->global_params);
		}
		//ApplyState(sys,N,GlobalVars,analysis->ve,t_b,workspace);


		//Write the analysis to disk	//!!!! Switch with set analysis? !!!!
		for(i=0;i<my_N;i++)		//Output the initial values to the file
		{
			current = sys[my_sys[i]];

			if(current->print_time > 0.0)
			{
//printf("Writing: %u\n",current->ID);
//Print_Vector(current->list->head->y_approx);
				(current->disk_iterations)++;

				fsetpos(outputfile,&(current->pos));
				fwrite(&(current->next_save),sizeof(double),1,outputfile);
				for(j=0;j<GlobalVars->num_print;j++)
					fwrite(&(current->list->head->y_approx->ve[GlobalVars->print_indices[j]]),sizeof(double),1,outputfile);
				fgetpos(outputfile,&(current->pos));
				current->next_save += current->print_time;
			}
		}

		//!!!! M and B may be calculated wrong. But I think this is correct... !!!!
/*
		//Propagate B
		mm_mlt(K,H,temp_NN);
		for(i=0;i<all_states;i++)	//Compute temp_NN = I - KH
		{
			for(j=0;j<i;j++)		temp_NN->me[i][j] *= -1.0;
			temp_NN->me[i][i] = 1.0 - temp_NN->me[i][i];
			for(j=i+1;j<all_states;j++)	temp_NN->me[i][j] *= -1.0;
		}
		mm_mlt(temp_NN,B,A);
		mm_mlt(M,A,temp_NN);
		mmT_mlt(temp_NN,M,B);

if( (k == stophere || stophere == -1) )
{
printf("A is\n");
Print_Matrix(A);

//printf("Previous M is %f %u\n",t_b,k);
//Print_Matrix(M);
}
*/


		//Reset variational equation start
		for(i=0;i<N;i++)
			ReadInitData(GlobalVars->global_params,sys[i]->params,sys[i]->iparams,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type);

		//Compute next background
		//AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,1,outputfile);
		ExactSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,1,outputfile);

/*
if(k == stophere)
{
printf("\nOut k = %u t = %f\n",k,t_b);
for(i=0;i<N;i++)
{
	printf("ID = %u\n",sys[i]->ID);
	Print_Vector(sys[i]->list->head->y_approx);
}
//getchar();
}
*/

		//Set the next matrix M
/*
		for(i=0;i<M->m;i++)
		{
			for(j=0;j<i;j++)	M->me[i][j] = 0.0;
			M->me[i][j] = 1.0;
			for(j=i+1;j<M->n;j++)	M->me[i][j] = 0.0;
		}
*/

		//New
		for(i=0;i<N;i++)
		{
			current = sys[i];

			//From my link
			for(j=0;j<(problem_dim-1);j++)
				M->me[i*problem_dim + 1 + j][i*problem_dim + 1 + j] = current->list->tail->y_approx->ve[problem_dim + j];

			for(j=0;j<problem_dim;j++)
				M->me[i*problem_dim][i*problem_dim + j] = current->list->tail->y_approx->ve[problem_dim + (problem_dim-1) + j];

			//From parents
			unsigned int counter = 3*problem_dim - 1;
			for(j=0;j<current->numparents;j++)
			{
				for(l=0;l<current->numupstream[j];l++)
				{
					for(m=0;m<problem_dim;m++)
					{
						M->me[i*problem_dim][current->upstream[j][l]*problem_dim + m] = current->list->tail->y_approx->ve[counter];
						counter++;
					}
				}
			}
		}


/*
		//Old
		for(i=0;i<all_states;i++)
		{
			unsigned int adjustment = (i % problem_dim) * all_states;
			for(j=0;j<all_states;j++)
			{
				M->me[i][j] = sys[i/problem_dim]->list->tail->y_approx->ve[problem_dim + j + adjustment];
			}
		}
*/


if( (k == stophere || stophere == -1) )
{
printf("M is\n");
Print_Matrix(M);
//printf("%u is\n",sys[0]->ID);
//Print_Vector(sys[0]->list->tail->y_approx);
//printf("%u is\n",sys[1]->ID);
//Print_Vector(sys[1]->list->tail->y_approx);
//printf("%u is\n",sys[2]->ID);
//Print_Vector(sys[2]->list->tail->y_approx);

//getchar();
}


		//Propagate B
		mm_mlt(K,H,temp_NN);
		for(i=0;i<all_states;i++)	//Compute temp_NN = I - KH
		{
			for(j=0;j<i;j++)		temp_NN->me[i][j] *= -1.0;
			temp_NN->me[i][i] = 1.0 - temp_NN->me[i][i];
			for(j=i+1;j<all_states;j++)	temp_NN->me[i][j] *= -1.0;
		}
		mm_mlt(temp_NN,B,A);
		mm_mlt(M,A,temp_NN);
		mmT_mlt(temp_NN,M,B);

/*
for(i=0;i<all_states;i++)
B->me[i][i] += 1e-8;
*/

if( (k == stophere || stophere == -1) )
{
printf("A is\n");
Print_Matrix(A);
}

		//Cleanup sys
		for(i=0;i<N;i++)	//Set time to 0.0
		{
			current = sys[i];
			if(current->list != NULL)
			{
//printf("Id = %u\n",current->ID);
//Print_Vector(current->list->head->y_approx);
				while(current->current_iterations > 1)
				{
					Remove_Head_Node(current->list);
					(current->current_iterations)--;
				}
//Print_Vector(current->list->head->y_approx);
				current->list->head->t = GlobalVars->maxtime;
				current->last_t = GlobalVars->maxtime;
				current->steps_on_diff_proc = 1;
				current->iters_removed = 0;
				current->rejected = 0;
				if(current->numparents == 0)	current->ready = 1;
				else				current->ready = 0;
			}
		}

		GlobalVars->maxtime += data_inc;
		t_b += data_inc;
	}

	//Free some memory
	TransData_Free(my_data);
	m_free(H);
	m_free(B);
	m_free(R);
	v_free(d);
	v_free(analysis);
	m_free(K_temp);
	m_free(K);
	m_free(M);
	m_free(A);
	m_free(BHT);
	m_free(temp_dd1);
	m_free(temp_dd2);
	m_free(temp_NN);
	v_free(y_b);
	v_free(true_init);
	free(ipiv);
	free(y_0);
	for(i=0;i<numlinks;i++)
	{
		for(j=0;j<numsteps[i];j++)
			free(truesolution[i][j]);
		free(truesolution[i]);
	}
	free(truesolution);
	free(ids);
	free(numsteps);
	v_free(init_stddev);

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
	total_time += Process_Data(sys,GlobalVars,N,save_list,save_size,my_save_size,id_to_loc,assignments,NULL);

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


void ExactSolver(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,unsigned int my_max_nodes,UnivVars* GlobalVars,int* assignments,unsigned int** id_to_loc,TempStorage* workspace,ConnData* conninfo,TransData* my_data,short int print_flag,FILE* outputfile)
{
	unsigned int i,j;
	double t,precip;
	double t_0 = sys[0]->last_t;
	Link* link_i;
	VEC* y_0_0 = sys[0]->list->head->y_approx;
	VEC* y_0_1 = sys[1]->list->head->y_approx;
	VEC* y_0_2 = sys[2]->list->head->y_approx;
	VEC* ans_0 = sys[0]->list->head->next->y_approx;
	VEC* ans_1 = sys[1]->list->head->next->y_approx;
	VEC* ans_2 = sys[2]->list->head->next->y_approx;

	for(i=0;i<N;i++)
	{
		sys[i]->list->tail = sys[i]->list->head->next;
		sys[i]->current_iterations = 2;
	}

	//Compute solution
	for(t=t_0 + GlobalVars->print_time;t*0.999<=GlobalVars->maxtime;t+=GlobalVars->print_time)
	{
		for(i=0;i<sys[0]->rain->n_times;i++)
		{
			if(sys[0]->rain->rainfall[i][0] <= t && t < sys[0]->rain->rainfall[i+1][0])
			{
				precip = sys[0]->rain->rainfall[i][1];
				break;
			}	
		}

		ExactYConn(t,t_0,y_0_0,y_0_1,y_0_2,sys[0]->params,sys[1]->params,sys[2]->params,GlobalVars->global_params,precip,ans_0,ans_1,ans_2);

		for(i=0;i<N;i++)
		{
			link_i = sys[i];
			link_i->last_t = t;

			if(print_flag)
			{
//if(i == 0)
//printf("%e %e %e %e\n",t,link_i->next_save,link_i->last_t,link_i->print_time);
				//while(t <= link_i->next_save && link_i->next_save <= link_i->last_t)
				if(t < GlobalVars->maxtime && t*0.999 <= link_i->next_save && link_i->next_save <= t*1.001)
				//if(t < GlobalVars->maxtime && fabs(t - link_i->next_save) < 1e-12)
				{
					//Don't write anything if using data assimilation and at a time when data is available
					if(fabs(GlobalVars->maxtime - link_i->next_save) < 1e-12)	break;
					//double rounded = 1e-13*rint(1e13*(GlobalVars->maxtime - link_i->next_save));
					//if(rounded < 1e-13 && -rounded < 1e-13)		break;

					(link_i->disk_iterations)++;

					//Write to a file
					fsetpos(outputfile,&(link_i->pos));
					fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
					for(j=0;j<GlobalVars->num_print;j++)
						fwrite(&(link_i->list->head->next->y_approx->ve[GlobalVars->print_indices[j]]),sizeof(double),1,outputfile);
					fgetpos(outputfile,&(link_i->pos));
					link_i->next_save += link_i->print_time;
				}
			}
		}
	}
}

//Solution at all links, for a y connection
//Works with model 15 and 315
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
//The numbering is:        0      1        2     3   4   5
void ExactYConn(double t,double t_0,VEC* y_0_0,VEC* y_0_1,VEC* y_0_2,VEC* params_0,VEC* params_1,VEC* params_2,VEC* global_params,double precip,VEC* ans_0,VEC* ans_1,VEC* ans_2)
{
	double tau_0 = 1.0/params_0->ve[12];
	double tau_1 = 1.0/params_1->ve[12];
	double tau_2 = 1.0/params_2->ve[12];

	double alpha_0 = params_0->ve[14];
	double alpha_1 = params_1->ve[14];
	double alpha_2 = params_2->ve[14];

	double beta_0 = params_0->ve[16];
	double beta_1 = params_1->ve[16];
	double beta_2 = params_2->ve[16];

	double gamma_0 = params_0->ve[17];
	double gamma_1 = params_1->ve[17];
	double gamma_2 = params_2->ve[17];

	double A_0 = (beta_0 * precip * tau_0 * alpha_0 - y_0_0->ve[1] * alpha_0) / (1.0 - tau_0 * gamma_0) + y_0_0->ve[0];
	double A_1 = (beta_1 * precip * tau_1 * alpha_1 - y_0_1->ve[1] * alpha_1) / (1.0 - tau_1 * gamma_1) + y_0_1->ve[0];	
	double B_0 = (y_0_0->ve[1] * gamma_0 * alpha_0 - beta_0 * precip * alpha_0) / ((1.0 - tau_0 * gamma_0) * gamma_0);
	double B_1 = (y_0_1->ve[1] * gamma_1 * alpha_1 - beta_1 * precip * alpha_1) / ((1.0 - tau_1 * gamma_1) * gamma_1);
	double C_0 = beta_0 * precip * alpha_0 / gamma_0;
	double C_1 = beta_1 * precip * alpha_1 / gamma_1;

	double time_diff = t-t_0;
	double exp_tau_0 = exp(-1.0/tau_0 * time_diff);
	double exp_tau_1 = exp(-1.0/tau_1 * time_diff);
	double exp_tau_2 = exp(-1.0/tau_2 * time_diff);
	double exp_gamma_0 = exp(-gamma_0 * time_diff);
	double exp_gamma_1 = exp(-gamma_1 * time_diff);
	double exp_gamma_2 = exp(-gamma_2 * time_diff);

	//Solution to model eqs
	ans_0->ve[0] = A_0 * exp_tau_0 + B_0 * exp_gamma_0 + C_0;
	ans_0->ve[1] = beta_0 * precip / gamma_0 - (beta_0 * precip / gamma_0 - y_0_0->ve[1]) * exp_gamma_0;
	ans_1->ve[0] = A_1 * exp_tau_1 + B_1 * exp_gamma_1 + C_1;
	ans_1->ve[1] = beta_1 * precip / gamma_1 - (beta_1 * precip / gamma_1 - y_0_1->ve[1]) * exp_gamma_1;

	if( fabs(tau_0 - tau_2) < 1e-12 && fabs(tau_1 - tau_2) < 1e-12 )	//!!!! Ok, there are other possibilities... !!!!
	{	//The A_0 and A_1 terms are treated differently here
		ans_2->ve[0] = beta_2 * precip * alpha_2 / gamma_2 
			- (alpha_2 * beta_2 * precip - y_0_2->ve[1] * gamma_2 * alpha_2)/((1.0 - tau_2 * gamma_2) * gamma_2) * exp_gamma_2
			+ 1.0/tau_0 * A_0 * exp_tau_0 * time_diff + B_0/(1.0-tau_2*gamma_0)*exp_gamma_0 + C_0
			+ 1.0/tau_1 * A_1 * exp_tau_1 * time_diff + B_1/(1.0-tau_2*gamma_1)*exp_gamma_1 + C_1
			+ y_0_2->ve[0] * exp_tau_2
			- ( beta_2 * precip * alpha_2 / gamma_2 - (alpha_2 * beta_2 * precip - y_0_2->ve[1] * gamma_2 * alpha_2)/((1.0-tau_2*gamma_2) * gamma_2)
				+ B_0/(1.0-tau_2*gamma_0) + C_0
				+ B_1/(1.0-tau_2*gamma_1) + C_1 ) * exp_tau_2;
	}
	else
	{
		ans_2->ve[0] = beta_2 * precip * alpha_2 / gamma_2 
			- (alpha_2 * beta_2 * precip - y_0_2->ve[1] * gamma_2 * alpha_2)/((1.0 - tau_2 * gamma_2) * gamma_2) * exp_gamma_2
			+ A_0/(1.0-tau_2/tau_0) * exp_tau_0 + B_0/(1.0-tau_2*gamma_0)*exp_gamma_0 + C_0
			+ A_1/(1.0-tau_2/tau_1) * exp_tau_1 + B_1/(1.0-tau_2*gamma_1)*exp_gamma_1 + C_1
			+ y_0_2->ve[0] * exp_tau_2
			- ( beta_2 * precip * alpha_2 / gamma_2 - (alpha_2 * beta_2 * precip - y_0_2->ve[1] * gamma_2 * alpha_2)/((1.0-tau_2*gamma_2) * gamma_2)
				+ A_0/(1.0-tau_2/tau_0) + B_0/(1.0-tau_2*gamma_0) + C_0
				+ A_1/(1.0-tau_2/tau_1) + B_1/(1.0-tau_2*gamma_1) + C_1 ) * exp_tau_2;
	}

	ans_2->ve[1] = beta_2 * precip / gamma_2 - (beta_2 * precip / gamma_2 - y_0_2->ve[1]) * exp_gamma_2;

	//Solution to variational eqs
	ans_0->ve[2] = exp_gamma_0;
	ans_0->ve[3] = exp_tau_0;
	ans_0->ve[4] = alpha_0 / (1.0 - tau_0 * gamma_0) * (exp_gamma_0 - exp_tau_0);

	ans_1->ve[2] = exp_gamma_1;
	ans_1->ve[3] = exp_tau_1;
	ans_1->ve[4] = alpha_1 / (1.0 - tau_1 * gamma_1) * (exp_gamma_1 - exp_tau_1);

	if( fabs(tau_0 - tau_2) < 1e-12 && fabs(tau_1 - tau_2) < 1e-12 )	//!!!! Ok, there are other possibilities... !!!!
	{
		ans_2->ve[2] = exp_gamma_2;	//ds_2/ds_2i
		ans_2->ve[3] = exp_tau_2;	//dq_2/dq_2i
		ans_2->ve[4] = alpha_2 / (1.0 - tau_2 * gamma_2) * (exp_gamma_2 - exp_tau_2);	//dq_2/ds_2i
		ans_2->ve[5] = 1.0/tau_0 * exp_tau_0 * time_diff;	//dq_2/dq_0i
		ans_2->ve[6] = -alpha_0/(tau_0*(1.0-tau_0*gamma_0)) * exp_tau_0 * time_diff + alpha_0/((1.0-tau_0*gamma_0)*(1.0-tau_2*gamma_0))*(exp_gamma_0 - exp_tau_2);	//dq_2/ds_0i
		ans_2->ve[7] = 1.0/tau_0 * exp_tau_1 * time_diff;	//dq_2/dq_1i
		ans_2->ve[8] = -alpha_1/(tau_1*(1.0-tau_1*gamma_1)) * exp_tau_1 * time_diff + alpha_1/((1.0-tau_1*gamma_1)*(1.0-tau_2*gamma_1))*(exp_gamma_1 - exp_tau_2);	//dq_2/ds_1i
	}
	else
	{
		ans_2->ve[2] = exp_gamma_2;	//ds_2/ds_2i
		ans_2->ve[3] = exp_tau_2;	//dq_2/dq_2i
		ans_2->ve[4] = alpha_2 / (1.0 - tau_2 * gamma_2) * (exp_gamma_2 - exp_tau_2);	//dq_2/ds_2i
		ans_2->ve[5] = 1.0 / (1.0 - tau_2/tau_0) * (exp_tau_0 - exp_tau_2);	//dq_2/dq_0i
		ans_2->ve[6] = -alpha_0/((1.0-tau_0*gamma_0)*(1.0-tau_2/tau_0))*(exp_tau_0 - exp_tau_2) + alpha_0/((1.0-tau_0*gamma_0)*(1.0-tau_2*gamma_0))*(exp_gamma_0 - exp_tau_2);	//dq_2/ds_0i
		ans_2->ve[7] = 1.0 / (1.0 - tau_2/tau_1) * (exp_tau_1 - exp_tau_2);	//dq_2/dq_1i
		ans_2->ve[8] = -alpha_1/((1.0-tau_1*gamma_1)*(1.0-tau_2/tau_1))*(exp_tau_1 - exp_tau_2) + alpha_1/((1.0-tau_1*gamma_1)*(1.0-tau_2*gamma_1))*(exp_gamma_1 - exp_tau_2);	//dq_2/ds_1i
	}

/*
if(t > 60.0)
{
printf("********* %e\n",precip);
printf("%e\n",time_diff);
Print_Vector(ans_2);
getchar();
}
*/

}



double FindDischargeExact(double t,Link** sys,UnivVars* GlobalVars)
{
	double t_0 = 0.0;
	double precip = 0.0;
	VEC* y_0_0 = v_get(sys[0]->list->head->y_approx->dim);
	VEC* y_0_1 = v_get(sys[1]->list->head->y_approx->dim);
	VEC* y_0_2 = v_get(sys[2]->list->head->y_approx->dim);
	VEC* ans_0 = sys[0]->list->head->next->y_approx;
	VEC* ans_1 = sys[1]->list->head->next->y_approx;
	VEC* ans_2 = sys[2]->list->head->next->y_approx;

	y_0_0->ve[0] = y_0_1->ve[0] = y_0_2->ve[0] = 10.0;
	y_0_0->ve[1] = y_0_1->ve[1] = y_0_2->ve[1] = 0.0;

	ExactYConn(t,t_0,y_0_0,y_0_1,y_0_2,sys[0]->params,sys[1]->params,sys[2]->params,GlobalVars->global_params,precip,ans_0,ans_1,ans_2);

	v_free(y_0_0);
	v_free(y_0_1);
	v_free(y_0_2);

	return ans_2->ve[0];
}



/*
double ExactType0(double t,VEC* y_0,VEC* params,VEC* global_params,double stddev)
{
	double lambda_1 = global_params->ve[1];
	double invtau = params->ve[12];
	double TWO_PI = 2.0 * 3.141592653589;
	double mean = pow(1.0 + lambda_1 * invtau * pow(y_0->ve[0],lambda_1) * t,-1.0/lambda_1) * y_0->ve[0];

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
*/

