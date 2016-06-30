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
	VEC* init_stddev = v_get(problem_dim);		//The standard deviation of the init conditions
	init_stddev->ve[0] = 1e-1;
	init_stddev->ve[1] = 1e-1;
	//init_stddev->ve[0] = 1e+4;
	//init_stddev->ve[1] = 1e+3;
	double data_inc = 5.0;		//Time increment when new data becomes available
	//unsigned int data_locs[] = {0};	//Locations of links with data	!!!! Make sure these are locations and not IDs !!!!
	//unsigned int data_locs[] = {1,3,4};
	unsigned int data_locs[] = {2};
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

	//Record what entries of M are not always 0
	//if(GlobalVars->type != 0 && GlobalVars->type != 1)	printf("Warning in assimExtendedApproxSerial.c. Cannot set Mnot0 for type %u.\n",GlobalVars->type);
	unsigned int* Mnot0_size = (unsigned int*) malloc(all_states*sizeof(unsigned int));
	unsigned int idx,loc;
	for(i=0;i<all_states;i++)	Mnot0_size[i] = 1;
	for(i=0;i<N;i++)
	{
		current = sys[i];
		idx = i*problem_dim;
		for(j=0;j<current->numparents;j++)
			Mnot0_size[idx] += current->numupstream[j];
		Mnot0_size[idx] *= problem_dim;
	}

	unsigned int** Mnot0 = (unsigned int**) malloc(all_states*sizeof(unsigned int*));
	for(i=0;i<all_states;i++)	Mnot0[i] = (unsigned int*) malloc(Mnot0_size[i] * sizeof(unsigned int));

	for(i=0;i<all_states;i+=problem_dim)
	{
		l = 0;
		loc = i / problem_dim;
		for(j=0;j<sys[loc]->numparents;j++)
		{
			for(k=0;k<sys[loc]->numupstream[j];k++)
			{
				for(m=0;m<problem_dim;m++)
					Mnot0[i][l+m] = sys[loc]->upstream[j][k] * problem_dim + m;
				l += m;
			}
		}
		Mnot0[i][l] = loc*problem_dim;
		if(problem_dim > 1)	Mnot0[i][l+1] = loc*problem_dim + 1;
		l += problem_dim;
		if(problem_dim > 1)	Mnot0[i+1][0] = i+1;
		merge_sort_1D(Mnot0[i],l);	//!!!! Not sure if this is needed !!!!
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
		for(l=0;l<problem_dim;l++)	B->me[i+l][i+l] = sq(init_stddev->ve[l]);
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
	array[2*i+1] = 0.0;
}
	
for(i=0;i<N;i++)
	for(j=0;j<problem_dim;j++)
		sys[i]->list->head->y_approx->ve[j] = array[i*problem_dim+j];


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

/*
for(i=0;i<N;i++)
{
	printf("%i\n",i);
	Print_Vector(sys[i]->params);
	printf("*******\n");
}
getchar();
*/

double stophere = -60.0;
flaggy = -1;


	//Main loop
	for(k=0;k<iterations;k++)
	{
		t_b += data_inc;

printf("\n******************************************\n");
printf("t_b = %f k = %u / %u  Maxtime = %f at %f\n",t_b,k,iterations,GlobalVars->maxtime,sys[0]->last_t);


		//Calculate A
		mm_mlt(M,B,temp_NN);
		mmT_mlt(temp_NN,M,A);

if(t_b > stophere)
{
printf("A is\n");
Print_Matrix(A);
}

/*
if(t_b > stophere)
{
printf("R is\n");
Print_Matrix(R);
printf("H is \n");
Print_Matrix(H);
}
*/
		//Calculate gain matrix
		mmT_mlt(A,H,BHT);			// !!!! Simplify because of H? !!!!
		mm_mlt(H,BHT,temp_dd2);		// !!!! Simplify because of H? !!!!
		for(i=0;i<R->m;i++)	temp_dd2->me[i][i] += R->me[i][i];
		clapack_dgetrf(CblasRowMajor,numdata,numdata,temp_dd2->array,temp_dd2->m,ipiv);
		clapack_dgetri(CblasRowMajor,numdata,temp_dd2->array,temp_dd2->m,ipiv);
		mm_mlt(BHT,temp_dd2,K);		//!!!! Need to use ipiv? !!!!



		//Reset variational equation start
		for(i=0;i<N;i++)
			ReadInitData(GlobalVars->global_params,sys[i]->params,sys[i]->iparams,NULL,0,sys[i]->list->head->y_approx,GlobalVars->type);


/*
if(t_b > stophere)
{
printf("Going in at time %f\n",sys[0]->last_t);
for(i=0;i<N;i++)
{
printf("ID = %u\n",sys[i]->ID);
Print_Vector(sys[i]->list->tail->y_approx);
}
printf("***\n");
}
*/
		//Compute next background
		AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,1,outputfile);
/*
if(t_b > stophere)
{
printf("Out\n");
for(i=0;i<N;i++)
{
printf("ID = %u\n",sys[i]->ID);
Print_Vector(sys[i]->list->tail->y_approx);
}
printf("***\n");
}
*/
		//Set background
		for(i=0;i<N;i++)
			for(j=0;j<problem_dim;j++)
				y_b->ve[i*problem_dim+j] = sys[i]->list->tail->y_approx->ve[j];
				//y_b->ve[i*problem_dim+j] = sys[i]->list->head->y_approx->ve[j];

if(t_b > stophere)
{
printf("Background\n");
Print_Vector(y_b);
}


		//Calculate innovations
		for(i=0;i<numdata;i++)
		{
			current = sys[data_locs[i]];
			d->ve[i] = FindDischarge(truesolution,current->ID,t_b,numlinks,ids,numsteps,data_stddev) - current->list->tail->y_approx->ve[0];
/*
if(t_b > stophere)
{
printf("ID = %u  d = %e tail = %f  Got = %f\n",current->ID,d->ve[i],current->list->tail->y_approx->ve[0],FindDischarge(truesolution,current->ID,t_b,numlinks,ids,numsteps,data_stddev));
}
*/
		}


//for(i=0;i<numdata;i++)
//	d->ve[i] = 0.0;


if(t_b > stophere)
{
printf("Innovations\n");
Print_Vector(d);
}


		//Calculate analysis
		mv_mlt(K,d,analysis);
		v_add(y_b,analysis,analysis,0);		//!!!! Can this be rewritten to solve a linear system? !!!!


if(t_b > stophere)
{
printf("K is\n");
Print_Matrix(K);
printf("Analysis\n");
Print_Vector(analysis);
//getchar();
}

/*
if(t_b > stophere)
{
printf("Energy is %e\n",Functional(A,R,d,analysis,y_b));
}
*/

/*
if(t_b > stophere)
{
double storage = 0.0;
for(i=0;i<N;i++)
	storage += 60.0/sys[i]->params->ve[12] * y_b->ve[2*i] + sys[i]->params->ve[1] * y_b->ve[2*i+1];
printf("storage background = %f\n",storage);

storage = 0.0;
for(i=0;i<N;i++)
	storage += 60.0/sys[i]->params->ve[12] * analysis->ve[2*i] + sys[i]->params->ve[1] * analysis->ve[2*i+1];
printf("storage analysis = %f\n",storage);
}
*/


		//Propagate B
		mm_mlt(K,H,temp_NN);
		for(i=0;i<all_states;i++)	//Compute temp_NN = I - KH
		{
			for(j=0;j<i;j++)		temp_NN->me[i][j] *= -1.0;
			temp_NN->me[i][i] = 1.0 - temp_NN->me[i][i];
			for(j=i+1;j<all_states;j++)	temp_NN->me[i][j] *= -1.0;
		}
		mm_mlt(temp_NN,A,B);

/*
for(i=0;i<N;i++)
{
	B->me[i*2][i*2] += (B->me[i*2][i*2] >= 0.0) ? 1e+2 : -1e+2;
	B->me[i*2+1][i*2+1] += (B->me[i*2+1][i*2+1] >= 0.0) ? 1e+4 : -1e+4;
	//B->me[i*2][i*2] *= 1e4;
	//B->me[i*2+1][i*2+1] *= 1e4;
	//B->me[i*2][i*2] += 1e-2;
	//B->me[i*2+1][i*2+1] += 1e-2;
}
*/


if(t_b > stophere)
{
printf("B is\n");
Print_Matrix(B);
//getchar();
}




		//Set the next matrix M

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

/*
if(t_b > stophere)
{
printf("Before\n");
Print_Matrix(M);
}
*/

/*
//Adjust M
for(i=0;i<all_states;i++)
{
	for(j=0;j<Mnot0_size[i];j++)
		if(M->me[i][Mnot0[i][j]] < 1e-14)	M->me[i][Mnot0[i][j]] = 1e-14;
}
*/

/*
printf("After\n");
Print_Matrix(M);
getchar();
*/


if(t_b > stophere)
{
printf("M is %f\n",t_b);
Print_Matrix(M);
//getchar();
}



		//Reset sys
		for(i=0;i<N;i++)	//Set time to 0.0
		{
			current = sys[i];
			if(current->list != NULL)
			{
				while(current->current_iterations > 1)
				{
					Remove_Head_Node(current->list);
					(current->current_iterations)--;
				}
				current->list->head->t = GlobalVars->maxtime;
				current->last_t = GlobalVars->maxtime;
				current->steps_on_diff_proc = 1;
				current->iters_removed = 0;
				current->rejected = 0;
				if(current->numparents == 0)	current->ready = 1;
				else				current->ready = 0;
			}
		}


		//Set analysis
		for(i=0;i<N;i++)
		{
			for(j=0;j<problem_dim;j++)
				sys[i]->list->head->y_approx->ve[j] = analysis->ve[i*problem_dim+j];
			CheckConsistency(sys[i]->list->head->y_approx,GlobalVars->dim,GlobalVars->type,sys[i]->params,GlobalVars->global_params);
		}

		//Write the analysis to disk	//!!!! Switch with set analysis? !!!!
		for(i=0;i<my_N;i++)		//Output the initial values to the file
		{
			current = sys[my_sys[i]];
			if(current->print_time > 0.0 && (current->next_save * 0.999999 < current->last_t && current->last_t < current->next_save * 1.000001))
			{
//if(current->ID == 0)
//printf("Writing at %f  last_t = %f\n",current->next_save,current->last_t);
				(current->disk_iterations)++;

				fsetpos(outputfile,&(current->pos));
				fwrite(&(current->next_save),sizeof(double),1,outputfile);
				for(j=0;j<GlobalVars->num_print;j++)
					fwrite(&(current->list->head->y_approx->ve[GlobalVars->print_indices[j]]),sizeof(double),1,outputfile);
				fgetpos(outputfile,&(current->pos));
				current->next_save += current->print_time;
			}
		}

		//Go to the next time
		GlobalVars->maxtime += data_inc;
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
	for(i=0;i<all_states;i++)
		free(Mnot0[i]);
	free(Mnot0);
	free(Mnot0_size);

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

