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
#include "mathmethods.h"

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
	unsigned int *my_sys,my_N,my_max_nodes,N,nummethods,my_save_size,save_size,i,j,k,l,m,*save_list,peaksave_size,loc,idx;
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

	if(GlobalVars->type >= 300)
	{
		if(my_rank == 0)	printf("Error: Type %u is not allowed in ExtendedApproxSerial.\n",GlobalVars->type);
		MPI_Finalize();
		return 1;
	}

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
	//unsigned int problem_dim = GlobalVars->problem_dim;
	unsigned int dim = GlobalVars->dim;
	double data_stddev = 1e-4;		//The standard deviation of the data
	VEC* init_stddev = v_get(dim);		//The standard deviation of the init conditions
	init_stddev->ve[0] = 1e+2;
	init_stddev->ve[1] = 1e+1;
	double data_inc = 5.0;		//Time increment when new data becomes available
	//VEC* true_init = v_get(GlobalVars->dim);
	//true_init->ve[0] = 10.0;
	double delta_x_nom = 1e-6;
	double delta_x;
	//unsigned int data_locs[] = {0,1};	//Locations of links with data	!!!! Make sure these are locations and not IDs !!!!
	//unsigned int data_locs[] = {1,3,4};
	unsigned int data_locs[] = {0};
	unsigned int numdata = sizeof(data_locs)/sizeof(unsigned int);	//Number of links for which data is available
	double t_f = GlobalVars->maxtime;	//The actual final time. This gets overwritten in GlobalVars
	GlobalVars->maxtime = data_inc;
	double t_b;			//The background time
	unsigned int all_states = N*dim;
	MAT* H = m_get(numdata,all_states);
	MAT* B = m_get(all_states,all_states);
	MAT* R = m_get(numdata,numdata);
	VEC* d = v_get(numdata);	//Innovations
	VEC* analysis = v_get(all_states);	//Memory to hold the analysis
	MAT* K_temp = m_get(numdata,numdata);
	MAT* K = m_get(all_states,numdata);
	MAT* MT = m_get(all_states,all_states);
	MAT* A = m_get(all_states,all_states);
	MAT* BHT = m_get(all_states,numdata);
	MAT* temp_dd1 = m_get(numdata,numdata);
	MAT* temp_dd2 = m_get(numdata,numdata);
	MAT* temp_NN = m_get(all_states,all_states);
	VEC* y_b = v_get(all_states);
	int* ipiv = (int*) malloc(numdata*sizeof(int));
	double* y_0 = (double*) malloc(all_states*sizeof(double));
	VEC* true_init = v_get(all_states);
	for(i=0;i<all_states;i+=dim)
	{
		if(assignments[i/dim] == my_rank)
		{
			for(j=0;j<dim;j++)
				true_init->ve[i+j] = sys[i/dim]->list->head->y_approx->ve[j];
		}
	}
	VEC* next_sys_state = v_get(all_states);
	VEC* prev_sys_state = v_get(all_states);
	VEC* dummyvec = NULL;

	//Barrier and stop clock
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time = difftime(stop,start);
	if(my_rank == 0)	printf("Finished reading files. Total time: %f\n",total_time);

	//Compute which entries of MT are not always 0
	if(GlobalVars->type != 0 && GlobalVars->type != 1)	printf("Warning in assimExtendedApproxSerial.c. Cannot set Mnot0 for type %u.\n",GlobalVars->type);
	unsigned int* Mnot0_size = (unsigned int*) malloc(all_states*sizeof(unsigned int));
	for(i=0;i<all_states;i++)	Mnot0_size[i] = 1;
	for(i=0;i<N;i++)
	{
		current = sys[i];
		idx = i*dim;
		for(j=0;j<current->numparents;j++)
			Mnot0_size[idx] += current->numupstream[j];
		Mnot0_size[idx] *= dim;
	}

	unsigned int** Mnot0 = (unsigned int**) malloc(all_states*sizeof(unsigned int*));
	for(i=0;i<all_states;i++)	Mnot0[i] = (unsigned int*) malloc(Mnot0_size[i] * sizeof(unsigned int));

	for(i=0;i<all_states;i+=dim)
	{
		l = 0;
		loc = i / dim;
		for(j=0;j<sys[loc]->numparents;j++)
		{
			for(k=0;k<sys[loc]->numupstream[j];k++)
			{
				for(m=0;m<dim;m++)
					Mnot0[i][l+m] = sys[loc]->upstream[j][k] * dim + m;
				l += m;
			}
		}
		Mnot0[i][l] = loc*dim;
		if(dim > 1)	Mnot0[i][l+1] = loc*dim + 1;
		l += dim;
		if(dim > 1)	Mnot0[i+1][0] = i+1;
		merge_sort_1D(Mnot0[i],l);	//!!!! Not sure if this is needed !!!!
	}
/*
printf("Sparsity of M is\n");
for(i=0;i<all_states;i++)
{
	for(j=0;j<Mnot0_size[i];j++)
		printf("%u ",Mnot0[i][j]);
	printf("\n");
}
getchar();
*/
	//Compute which entries of M are not always 0
	unsigned int* MTnot0_size = (unsigned int*) calloc(all_states,sizeof(unsigned int));
	for(i=0;i<all_states;i++)
		for(j=0;j<Mnot0_size[i];j++)	MTnot0_size[Mnot0[i][j]]++;
	unsigned int** MTnot0 = (unsigned int**) malloc(all_states*sizeof(unsigned int*));
	for(i=0;i<all_states;i++)	MTnot0[i] = (unsigned int*) malloc(MTnot0_size[i]*sizeof(unsigned int));

	unsigned int* counter = (unsigned int*) calloc(all_states,sizeof(unsigned int));
	unsigned int row;
	for(i=0;i<all_states;i++)
	{
		for(j=0;j<Mnot0_size[i];j++)
		{
			row = Mnot0[i][j];
			MTnot0[row][counter[row]] = i;
			counter[row]++;
		}
	}
/*
printf("Sparsity of MT is\n");
for(i=0;i<all_states;i++)
{
	for(j=0;j<MTnot0_size[i];j++)
		printf("%u ",MTnot0[i][j]);
	printf("\n");
}
getchar();
*/
	//Clean up
	for(i=0;i<all_states;i++)	free(Mnot0[i]);
	free(Mnot0);
	free(Mnot0_size);
	free(counter);

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time += difftime(stop,start);
	if(my_rank == 0)	printf("Finished building sparsity structure. Total time: %f\n",difftime(stop,start));


	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,my_N);
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);


	//Initialize matrices
	for(i=0;i<numdata;i++)	H->me[i][data_locs[i]*dim] = 1.0;
	for(i=0;i<all_states;i+=dim)
		for(l=0;l<dim;l++)	B->me[i+l][i+l] = sq(init_stddev->ve[l]);
	for(i=0;i<numdata;i++)	R->me[i][i] = sq(data_stddev);
	for(i=0;i<all_states;i++)	MT->me[i][i] = 1.0;

	//Set init conditions
	double TWO_PI = 2.0 * 3.141592653589;
	double randnum1,randnum2;
	//for(i=0;i<my_N;i++)
	//	sys[my_sys[i]]->list->head->y_approx->ve[0] = true_init->ve[0] + init_stddev * sqrt(-2.0*log( (rand()%100000)/100000.0 + .00000001 ))*cos(TWO_PI * ((rand()%100000)/100000.0 + .00000001));
	for(i=0;i<my_N*dim;i+=dim)
	{
		//if(assignments[i/dim] == my_rank)
		{
			for(l=0;l<dim;l++)
			{
				randnum1 = (rand()%100000)/100000.0 + .00000001;
				randnum2 = ((rand()%100000)/100000.0 + .00000001);
				sys[my_sys[i/dim]]->list->head->y_approx->ve[l] = max( true_init->ve[dim*sys[my_sys[i/dim]]->location + l] + init_stddev->ve[l] * sqrt(-2.0*log( randnum1 ))*cos(TWO_PI * randnum2) , 0.0 );
			}
		}
	}
	GetState(sys,my_sys,N,my_N,assignments,GlobalVars,prev_sys_state->ve);

/*
//Print_VectorC(prev_sys_state);
//double array[] = { 16.5275766899834906, 0.1250310549029123, 11.2495128939127813, 0.9984767446669658, 7.4678815629879995, 0.3631496938064627, 0.0000000000000000, 0.0000000000000000, 0.3524852009873456, 0.0000000000000000, 3.7464509476927832, 0.0000000000000000, 14.6253464504624873, 0.0000000000000000, 19.1263612558515348, 0.6549741437215924, 0.0000000000000000, 0.0000000000000000, 2.4287923338950161, 0.0000000000000000, 31.9581627447345227, 0.3046019826458240 };
//double array[] = { 0.3524852009873456, 0.0000000000000000, 19.1263612558515348, 0.6549741437215924 };
double array[] = { 21.3084396756748049, 0.0000000000000000, 0.9646386156869209, 0.0000000000000000, 11.4054449810795546, 0.2450490919781408, 17.9312794563208122, 0.0000000000000000, 13.8679425051719711, 0.4730965170219442, 4.7993763516032608, 0.9378456156897189, 15.1712554944539981, 0.3319546970447389, 17.9569645539304616, 0.0000000000000000, 4.8853368410571489, 0.0000000000000000, 27.7140292609099816, 1.9564235223044673, 21.9937744492977103, 0.0000000000000000 };
//double array[] = { 21.3084396756748049, 0.0000000000000000, 4.7993763516032608, 0.9378456156897189 };
*/
double array[] = { 10.0, 0.0, 10.0, 0.0, 10.0, 0.0 };
for(i=0;i<all_states;i++)
	prev_sys_state->ve[i] = array[i];
ApplyState(sys,N,GlobalVars,prev_sys_state->ve,0.0);


	//Read in the "true" solution
	unsigned int numlinks,*ids,*numsteps;
	double*** truesolution = ReadSolution("TempData/testobservations.dat",&numlinks,&ids,&numsteps);

	//Calculate number of times to loop
	unsigned int iterations = rint((t_f - GlobalVars->t_0) / data_inc);
	t_b = GlobalVars->t_0;

unsigned int stophere = -2;
flaggy = -1;

	//Main loop
	//for(t_b = GlobalVars->t_0; t_b <= t_f; t_b += data_inc)
	for(k=0;k<iterations;k++)
	{
		//Set background
		for(i=0;i<N;i++)
			for(j=0;j<dim;j++)
				y_b->ve[i*dim+j] = sys[i]->list->head->y_approx->ve[j];

if(k == stophere || stophere == 500)
{
printf("Background\n");
Print_Vector(y_b);
}

		//Calculate innovations
		for(i=0;i<numdata;i++)
		{
			current = sys[data_locs[i]];
			d->ve[i] = FindDischarge(truesolution,current->ID,t_b,numlinks,ids,numsteps,data_stddev) - current->list->tail->y_approx->ve[0];
			//d->ve[i] = ExactType0(t_b,true_init,current->params,GlobalVars->global_params,data_stddev) - current->list->tail->y_approx->ve[0];
			//d->ve[0] = ExactType0_yconn(t_b,true_init,true_init,true_init,current->params,current->parents[0]->params,current->parents[1]->params,GlobalVars->global_params,data_stddev) - current->list->tail->y_approx->ve[0];
			//d->ve[i] = ExactTestType0(data_locs[i],t_b,y_0,sys[my_sys[0]]->params->ve[12],GlobalVars->global_params,data_stddev) - current->list->tail->y_approx->ve[0];
/*
if(i == 0)
{
printf("Got %.16f, current is %.16f, diff is %.16f\n",FindDischarge(truesolution,current->ID,t_b,numlinks,ids,numsteps,data_stddev),current->list->tail->y_approx->ve[0],d->ve[i]);
getchar();
}
*/
		}

if(k == stophere || stophere == 500)
{
printf("Innovations\n");
Print_Vector(d);
}

		//Calculate gain matrix
/*
if(k == stophere || stophere == 500)
{
printf("B is\n");
Print_Matrix(B);
//printf("H is\n");
//Print_Matrix(H);
//printf("R is\n");
//Print_Matrix(R);
//getchar();
}
*/
		mmT_mlt(B,H,BHT);			// !!!! Simplify because of H? !!!!
		mm_mlt(H,BHT,temp_dd2);		// !!!! Simplify because of H? !!!!
		for(i=0;i<R->m;i++)	temp_dd2->me[i][i] += R->me[i][i];
		clapack_dgetrf(CblasRowMajor,numdata,numdata,temp_dd2->array,temp_dd2->m,ipiv);
		clapack_dgetri(CblasRowMajor,numdata,temp_dd2->array,temp_dd2->m,ipiv);
		mm_mlt(BHT,temp_dd2,K);		//!!!! Need to use ipiv? !!!!

if(k == stophere || stophere == 500)
//if(t_b > 60.0)
{
printf("K is %f\n",t_b);
Print_Matrix(K);
//getchar();
}

		//Calculate analysis
		mv_mlt(K,d,analysis);
		v_add(y_b,analysis,analysis,0);		//!!!! Can this be rewritten to solve a linear system? !!!!
if(k == stophere || stophere == 500)
{
printf("Analysis\n");
Print_Vector(analysis);
//getchar();
}

		//Set analysis
/*
		for(i=0;i<N;i++)
			for(j=0;j<dim;j++)
				sys[i]->list->head->y_approx->ve[j] = analysis->ve[i*dim+j];
*/
		ApplyState(sys,N,GlobalVars,analysis->ve,t_b);

		//Write the analysis to disk	//!!!! Switch with set analysis? !!!!
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
				fprintf(outputfile,"%i %.12e ",current->ID,current->last_t);
				for(j=0;j<GlobalVars->num_print;j++)
					fprintf(outputfile,"%.12e ",current->list->head->y_approx->ve[GlobalVars->dense_indices[j]]);
				fprintf(outputfile,"\n");
*/
				current->next_save += current->print_time;
			}
		}

		//Propagate B
		mm_mlt(K,H,temp_NN);
		for(i=0;i<all_states;i++)	temp_NN->me[i][i] = 1.0 - temp_NN->me[i][i];
		mm_mlt(temp_NN,B,A);
		//mm_mlt(MT,A,temp_NN);
		//mmT_mlt(temp_NN,MT,B);
		mTm_mlt(MT,A,temp_NN);
		mm_mlt(temp_NN,MT,B);

/*
if(k == stophere || stophere == 500)
{
printf("A is\n");
Print_Matrix(A);
getchar();
}
*/
/*
for(i=0;i<N;i++)
printf("ID = %u  %.16f\n",sys[i]->ID,sys[i]->list->head->y_approx->ve[1]);
getchar();
*/
		//Compute next background
		AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,1,outputfile);

		//Set the next matrix M
		//for(i=0;i<N;i++)
		//	for(j=0;j<N;j++)
		//		M->me[i][j] = sys[i]->list->tail->y_approx->ve[GlobalVars->problem_dim + j];

		GetState(sys,my_sys,N,my_N,assignments,GlobalVars,next_sys_state->ve);	//Backup current state
		for(i=0;i<all_states;i++)	//!!!! This is only working in serial, for 1 dim problem !!!!
		{
			loc = i / dim;
			idx = i % dim;

			ApplyState(sys,N,GlobalVars,prev_sys_state->ve,t_b);	//Apply the state at time t_b
			//delta_x = max(0.000004791 * sys[loc]->list->head->y_approx->ve[idx],delta_x_nom);
			delta_x = delta_x_nom;

			if(sys[loc]->list->head->y_approx->ve[idx] - delta_x >= 0.0)	//Use a central difference
			{
				sys[loc]->list->head->y_approx->ve[idx] += delta_x;	//Perturb one direction
				//CheckConsistency(sys[loc]->list->head->y_approx,dim,GlobalVars->type,sys[loc]->params,GlobalVars->global_params);
				AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,0,outputfile);	//Solve system !!!! Can do this better !!!!
				//for(j=0;j<all_states;j++)	MT->me[i][j] = sys[j/dim]->list->tail->y_approx->ve[j%dim];
				for(j=0;j<MTnot0_size[i];j++)	MT->me[i][MTnot0[i][j]] = sys[MTnot0[i][j]/dim]->list->tail->y_approx->ve[MTnot0[i][j]%dim];

				ApplyState(sys,N,GlobalVars,prev_sys_state->ve,t_b);	//Apply the state at time t_b
				sys[loc]->list->head->y_approx->ve[idx] -= delta_x;	//Perturb the other direction
				CheckConsistency(sys[loc]->list->head->y_approx,dim,GlobalVars->type,sys[loc]->params,GlobalVars->global_params);
				AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,0,outputfile);	//Solve system !!!! Can do this better !!!!
				//for(j=0;j<all_states;j++)	MT->me[i][j] = (MT->me[i][j] - sys[j/dim]->list->tail->y_approx->ve[j%dim])/(2.0*delta_x);
				for(j=0;j<MTnot0_size[i];j++)	MT->me[i][MTnot0[i][j]] = (MT->me[i][MTnot0[i][j]] - sys[MTnot0[i][j]/dim]->list->tail->y_approx->ve[MTnot0[i][j]%dim])/(2.0*delta_x);
			}
			else	//Use a forward difference
			{
				//delta_x *= 10.0;

				sys[loc]->list->head->y_approx->ve[idx] += delta_x;	//Perturb one direction
				//CheckConsistency(sys[loc]->list->head->y_approx,dim,GlobalVars->type,sys[loc]->params,GlobalVars->global_params);
				AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,0,outputfile);	//Solve system !!!! Can do this better !!!!
				//for(j=0;j<all_states;j++)	MT->me[i][j] = (sys[j/dim]->list->tail->y_approx->ve[j%dim] - next_sys_state->ve[(j/dim)*dim + (j%dim)])/delta_x;
				for(j=0;j<MTnot0_size[i];j++)	MT->me[i][MTnot0[i][j]] = (sys[MTnot0[i][j]/dim]->list->tail->y_approx->ve[MTnot0[i][j]%dim] - next_sys_state->ve[(MTnot0[i][j]/dim)*dim + (MTnot0[i][j]%dim)])/delta_x;
			}
		}
		ApplyState(sys,N,GlobalVars,next_sys_state->ve,t_b+data_inc);
		


if(k == stophere || stophere == 500)
//if(t_b > 60.0)
{
printf("MT is %f\n",t_b);
Print_Matrix(MT);
getchar();
}


		//Prepare for next iteration
		GlobalVars->maxtime += data_inc;
		t_b += data_inc;

		dummyvec = prev_sys_state;
		prev_sys_state = next_sys_state;
		next_sys_state = dummyvec;
		dummyvec = NULL;
	}

	//Free some memory
	TransData_Free(my_data);
	for(i=0;i<all_states;i++)	free(MTnot0[i]);
	free(MTnot0);
	free(MTnot0_size);
	m_free(H);
	m_free(B);
	m_free(R);
	v_free(d);
	v_free(analysis);
	m_free(K_temp);
	m_free(K);
	m_free(MT);
	m_free(A);
	m_free(BHT);
	m_free(temp_dd1);
	m_free(temp_dd2);
	m_free(temp_NN);
	v_free(y_b);
	v_free(true_init);
	free(ipiv);
	free(y_0);
	v_free(next_sys_state);
	v_free(prev_sys_state);
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

