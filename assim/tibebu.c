#include <stdio.h>
#include <unistd.h>
#include "mpi.h"
#include "asynch_interface.h"

int my_rank;
int np;

int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user);
void Set_Output_User_LinkID(asynchsolver* asynch);

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
	time_t start,stop;
	double total_time;
	asynchsolver* asynch;
	int i;

	start = time(NULL);

	if(my_rank == 0)
		printf("\nBeginning initialization...\n*****************************\n");
	MPI_Barrier(MPI_COMM_WORLD);

	//Init asynch object and the river network
	asynch = Asynch_Init(MPI_COMM_WORLD);
	Asynch_Parse_GBL(asynch,argv[1]);
	Asynch_Load_System(asynch);

	if(my_rank == 0)
	{
		printf("\nModel type is %u.\nGlobal parameters are:\n",asynch->GlobalVars->type);
		Print_Vector(asynch->GlobalVars->global_params);
		printf("\n");
	}

	//Setup output for link id, if needed
	int id_setup = Asynch_Check_Output(asynch,"LinkID");
	if(id_setup != -1)
	{
		Set_Output_User_LinkID(asynch);
		Asynch_Set_Output(asynch,"LinkID",ASYNCH_INT,(void (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) &Output_Linkid,NULL,0);
	}

	//Prepare output files
	Asynch_Prepare_Temp_Files(asynch);
	Asynch_Write_Current_Step(asynch);
	//Asynch_Prepare_Peakflow_Output(asynch);
	Asynch_Prepare_Output(asynch);

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,asynch->my_N);
	sleep(1);
	MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank == 0)
	{
		stop = time(NULL);
		total_time = difftime(stop,start);
		printf("Finished initialization. Total time: %f\n\n\nComputing solution at each link...\n************************************\n",total_time);
	}
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	//Perform the calculations
	int num_intervals = 3;
	int print_peaks[] = {1,0,1};
	double ends[] = {60.0, 120.0, 240.0};
	char peakflowname[1024],givenpeakflowname[1024];
	Asynch_Get_Peakflow_Output_Name(asynch,givenpeakflowname);
	if(my_rank == 0)	printf("Peakflow file from .gbl is %s.\n",givenpeakflowname);
	double original_sim_time = Asynch_Get_Total_Simulation_Time(asynch);

	time(&start);
	for(i=0;i<num_intervals;i++)
	{
		sprintf(peakflowname,"%s_%.0f",givenpeakflowname,ends[i]);
		if(print_peaks[i] && my_rank == 0)	printf("Calculating for peakflow file %s...\n",peakflowname);
		Asynch_Set_Peakflow_Output_Name(asynch,peakflowname);

		Asynch_Set_Total_Simulation_Time(asynch,ends[i]);
		Asynch_Set_Peakflow_Output_Name(asynch,peakflowname);
		if(print_peaks[i])	Asynch_Prepare_Peakflow_Output(asynch);
		Asynch_Reset_Peakflow_Data(asynch);
		Asynch_Advance(asynch,1);
		if(print_peaks[i])
		{
			if(Asynch_Create_Peakflows_Output(asynch) && my_rank == 0)	printf("Error creating peakflow output\n");
		}
	}

	Asynch_Set_Total_Simulation_Time(asynch,original_sim_time);
	Asynch_Advance(asynch,1);

	MPI_Barrier(MPI_COMM_WORLD);
	time(&stop);

	//Out information
	total_time += difftime(stop,start);
	if(my_rank == 0)	printf("\nComputations complete. Total time for calculations: %f\n",difftime(stop,start));
	if(asynch->sys[asynch->my_sys[0]]->c == NULL)
	{
		printf("[%i]: The solution at ID %i at time %.12f is\n",my_rank,asynch->sys[asynch->my_sys[0]]->ID,asynch->sys[asynch->my_sys[0]]->last_t);
		Print_Vector(asynch->sys[asynch->my_sys[0]]->list->tail->y_approx);
	}

	//Take a snapshot
	Asynch_Take_System_Snapshot(asynch,NULL);

	//Create output files
	Asynch_Create_Output(asynch,NULL);
	//Asynch_Create_Peakflows_Output(asynch);

	//Cleanup
	Asynch_Delete_Temporary_Files(asynch);
	Asynch_Free(asynch);
	return 0;
}



int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user)
{
	return ((Link*)user)->ID;
}

//!!!! Gross, but not sure how else to handle this. Maybe with a lot of interface functions? !!!!
void Set_Output_User_LinkID(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		sys[my_sys[i]]->output_user = (void*) sys[my_sys[i]];
}


