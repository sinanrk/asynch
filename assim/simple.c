#include "rkmethods.h"
#include "system.h"
#include "comm.h"
#include "riversys.h"
#include "processdata.h"
#include <time.h>
#include "mpi.h"
#include "simplemethods.h"

extern void solout_(void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*);

int my_rank;
int np;             //number of processes

unsigned int rejects;

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
	unsigned int *my_sys,my_N,my_max_nodes,N,nummethods,my_save_size,save_size,i,j,k,*save_list;
	int *assignments,ii;
	char filename[128],rkdfilename[128];
	FILE* outputfile = NULL;
	short int *getting = NULL;
	TransData* my_data = NULL;
	RKMethod** AllMethods;
	Link** sys;
	unsigned int** id_to_loc;
	ErrorData* GlobalErrors;
	UnivVars* GlobalVars;
	TempStorage* workspace;
	//PGconn* conn = NULL;
	ConnData* conninfo = CreateConnData();

	start = time(NULL);

	//Read in .gbl file
	GlobalVars = Read_Global_Data(argv[1],&GlobalErrors,conninfo,rkdfilename);
	if(GlobalVars == NULL)	return 1;

	// **** Use additional input arguments here ****
	printf("\n!!!! Warning: Overwritting parameters in .gbl file !!!!\n\n");
	strcpy(GlobalVars->identifier,argv[3]);
	GlobalVars->global_params->ve[6] = atof(argv[2]);
	// **** ****

	//Read in remaining data from files
	sys = Create_River_System_parallel(rkdfilename,&N,&my_sys,&my_N,&my_max_nodes,&my_data,&assignments,&getting,&AllMethods,&nummethods,GlobalVars,GlobalErrors,&save_list,&my_save_size,&save_size,&id_to_loc,conninfo,&workspace);
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

	sprintf(filename,"%s.dat",outputfilename);
	openoutputfile_(filename);
	printf("Writing output to file %s.\n",filename);

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time = difftime(stop,start);
	if(my_rank == 0)	printf("Finished reading files. Total time: %f\n",total_time);

	//Initialize remaining data
	int is_ok;
	short int* done = calloc(my_max_nodes,sizeof(short int));
	short int parentsval;
	unsigned int alldone;
	RKSolutionNode *roottail;
	Link* current = sys[my_sys[my_N-1]];
	unsigned int last_idx,curr_idx,around;
	unsigned int two_my_N = 2*my_N;
	short int rain_flag = (GlobalVars->type == 3 || GlobalVars->type == 5);

	//Initialize values for rainfall data
	unsigned int increment,first_file,last_file,passes;
	int maxfileindex;
	double maxtime,file_time;

	if(rain_flag)
	{
		increment = GlobalVars->increment;
		first_file = GlobalVars->first_file;
		last_file = GlobalVars->last_file;
		file_time = GlobalVars->file_time;
		passes = (last_file-first_file+1)/increment;
		if((last_file-first_file+1)%increment != 0)	passes++;
	}
	else
	{
		passes = 1;
		maxtime = GlobalVars->maxtime;
	}

	//Allocate space for calculations
	VEC* temp1 = v_get(GlobalVars->dim);
	VEC* temp2 = v_get(GlobalVars->dim);
	VEC*** temp_parent_approx = (VEC***) malloc(GlobalVars->max_s*sizeof(VEC**));
	for(i=0;i<GlobalVars->max_s;i++)
	{
		temp_parent_approx[i] = malloc(GlobalVars->max_parents*sizeof(VEC*));
		for(j=0;j<GlobalVars->max_parents;j++)
			temp_parent_approx[i][j] = v_get(GlobalVars->dim);
	}

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,my_N);
	MPI_Barrier(MPI_COMM_WORLD);

	unsigned int bigdim = N*GlobalVars->dim;
	unsigned int lwork = 8*bigdim+21+5*1;
	double* work = calloc(lwork,sizeof(double));
	unsigned int liwork = 21;
	int* iwork = calloc(liwork,sizeof(int));
	unsigned long int* ipar = malloc(5*sizeof(unsigned long int));
	ipar[0] = (unsigned long int) sys;
	ipar[1] = N;
	ipar[2] = (unsigned long int) GlobalVars;
	//ipar[3] = 0*2+1;
	//ipar[3] = 177260*2-1;
	ipar[3] = 2527*2+1;
	ipar[4] = 0;	//Number of rejections
	int idid;

	work[4] = -1.0;
	work[6] = .1;	//Initial step size
	iwork[3] = -1;
	//work[3] = 1000.0;

	//Output
	int iout = 2;
	iwork[4] = 1;
	iwork[20] = ipar[3];

	double t_0 = 0.0;
	double* y_0 = malloc(bigdim*sizeof(double));
	int itol = 0;
	double *abstols,*reltols;
	if(itol == 1)
	{
		abstols = malloc(bigdim*sizeof(double));
		reltols = malloc(bigdim*sizeof(double));
		for(i=0;i<bigdim;i+=2)
		{
			abstols[i] = GlobalErrors->abstol->ve[0];
			reltols[i] = GlobalErrors->reltol->ve[0];
		}
		for(i=1;i<bigdim;i+=2)
		{
			abstols[i] = GlobalErrors->abstol->ve[1];
			reltols[i] = GlobalErrors->reltol->ve[1];
		}		
	}

	printf("Using special init conditions...\n");
	for(i=1;i<bigdim;i+=2)	y_0[i] = 1.0 / N * GlobalVars->global_params->ve[9];
	//for(i=0;i<bigdim;i+=2)	y_0[i] = GlobalVars->global_params->ve[0] / sys[i/GlobalVars->dim]->params->ve[0] * y_0[i+1];
	for(i=0;i<bigdim;i+=2)	y_0[i] = (GlobalVars->global_params->ve[0] + GlobalVars->global_params->ve[6]) / sys[i/GlobalVars->dim]->params->ve[0] * y_0[i+1];

	//Start the clock
	start = time(NULL);

	if(GlobalVars->type == 3 || GlobalVars->type == 105)
	{
		if(GlobalVars->rain_flag == 2)
		{
			unsigned int increment,first_file,last_file,passes,k;
			int maxfileindex;
			double maxtime,file_time;

			increment = GlobalVars->increment;
			first_file = GlobalVars->first_file;
			last_file = GlobalVars->last_file;
			file_time = GlobalVars->file_time;
			passes = (last_file-first_file+1)/increment;
			if((last_file-first_file+1)%increment != 0)	passes++;

			for(k=0;k<passes;k++)
			{
				if(k == passes-1)	maxtime = GlobalVars->maxtime;
				else			maxtime = min(GlobalVars->maxtime,(k+1)*file_time*increment);
				maxfileindex = (int) min((double)first_file+(k+1)*increment,(double)last_file);

				Create_Rain_Data_Par(sys,N,my_N,GlobalVars,my_sys,assignments,GlobalVars->rain_filename,first_file+k*increment,maxfileindex,k*file_time*increment,file_time,id_to_loc,increment+1);

				if(itol == 0)
					dopri5_(&bigdim,&full_system1,&t_0,y_0,&maxtime,&(GlobalErrors->reltol->ve[0]),&(GlobalErrors->abstol->ve[0]),&itol,&solout_,&iout,work,&lwork,iwork,&liwork,NULL,ipar,&idid);
				else
					dopri5_(&bigdim,&full_system1,&t_0,y_0,&maxtime,reltols,abstols,&itol,&solout_,&iout,work,&lwork,iwork, &liwork,NULL,ipar,&idid);

				t_0 = maxtime;
				if(GlobalVars->maxtime <= t_0)	break;
			}
		}
		else if(GlobalVars->rain_flag == 4)
		{
			if(itol == 0)
				dopri5_(&bigdim,&full_system1,&t_0,y_0,&(GlobalVars->maxtime),&(GlobalErrors->reltol->ve[0]),&(GlobalErrors->abstol->ve[0]),&itol,&solout_,&iout,work,&lwork,iwork,&liwork,NULL,ipar,&idid);
			else
				dopri5_(&bigdim,&full_system1,&t_0,y_0,&(GlobalVars->maxtime),reltols,abstols,&itol,&solout_,&iout,work,&lwork,iwork, &liwork,NULL,ipar,&idid);
		}
	}
	else if(GlobalVars->type == 0)
		dopri5_(&bigdim,&full_system0,&t_0,y_0,&(GlobalVars->maxtime),&(GlobalErrors->reltol->ve[0]),&(GlobalErrors->abstol->ve[0]),&itol,&solout_,&iout,work,&lwork,iwork,&liwork,NULL,ipar,&idid);
	else if(GlobalVars->type == 2)
		dopri5_(&bigdim,&full_system2,&t_0,y_0,&(GlobalVars->maxtime),&(GlobalErrors->reltol->ve[0]),&(GlobalErrors->abstol->ve[0]),&itol,&solout_,&iout,work,&lwork,iwork,&liwork,NULL,ipar,&idid);
	else if(GlobalVars->type == 200)
		dopri5_(&bigdim,&full_system200,&t_0,y_0,&(GlobalVars->maxtime),&(GlobalErrors->reltol->ve[0]),&(GlobalErrors->abstol->ve[0]),&itol,&solout_,&iout,work,&lwork,iwork,&liwork,NULL,ipar,&idid);

	//Stop the clock
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time += difftime(stop,start);
	printf("Time for calculations: %f\n",difftime(stop,start));

	//Print some data
	printf("The solution at ID %u is %.16f. idid is %i\n",sys[0]->ID,y_0[0],idid);
	printf("Total rejections: %u\n",ipar[4]);

	//Cleanup
	ConnData_Free(conninfo);
	free(work);
	free(iwork);
	free(ipar);
	free(y_0);
	if(itol == 1)
	{
		free(reltols);
		free(abstols);
	}
	Destroy_ErrorData(GlobalErrors);
	v_free(temp1);
	v_free(temp2);
	for(i=0;i<GlobalVars->max_s;i++)
	{
		for(j=0;j<GlobalVars->max_parents;j++)	v_free(temp_parent_approx[i][j]);
		free(temp_parent_approx[i]);
	}
	free(temp_parent_approx);
	free(done);
	free(assignments);
	free(getting);
	for(i=0;i<N;i++)
		free(id_to_loc[i]);
	free(id_to_loc);
	closeoutputfile_();
	//if(my_save_size > 0)	fclose(outputfile);
	MPI_Barrier(MPI_COMM_WORLD);	//Very important: make sure all processes have flushed their buffers before data is processed.

	//Last bit of cleanup
	for(i=0;i<N;i++)	Destroy_Link(sys[i],GlobalVars->iter_limit,rkdfilename[0] != '\0',GlobalVars);
	//for(i=0;i<N;i++)	Destroy_Link(sys[i],GlobalVars->iter_limit,rkdfilename[0] != '\0');
	free(sys);
	free(my_sys);
	for(i=0;i<nummethods;i++)	Destroy_RKMethod(AllMethods[i]);
	free(AllMethods);
	free(GlobalVars);

	MPI_Finalize();
	return 0;
}

