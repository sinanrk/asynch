#include "asynch_interface.h"

int main(int argc,char* argv[])
{
	//Initialize MPI
	MPI_Init(&argc,&argv);	//Start MPI
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);	//Find this process's rank

	//Make sure enough a command line input is given. Otherwise, quit.
	if(argc < 5)
	{
		if(my_rank == 0)	printf("\nCommand line parameters required:\n\tA universal variable file (.gbl),\n\tfile (.ini/.uini/.rec) for initial condition,\n\ttimestamp (epoch) of first rainfall time,\n\ttimestamp (epoch) of last rainfall time.\n\n");
		MPI_Finalize();
		return 1;
	}

	//Initialize Asynch Solver
	char database_info[256];
	//sprintf(database_info,"dbname=rm_model host=s-iihr58.iihr.uiowa.edu user=automated_solver password=C5.pfest0");	//Use server 58 for topology data
	asynchsolver* asynch = Asynch_Init(MPI_COMM_WORLD);
	Asynch_Parse_GBL(asynch,argv[1]);
	Asynch_Set_Total_Simulation_Time(asynch,24.0*60.0);
	Asynch_Set_Init_File(asynch,argv[2]);
	Asynch_Load_System(asynch);

	//Set timestamps
	Asynch_Set_RainDB_Starttime(asynch,atoi(argv[3]),0);
	Asynch_Set_First_Rainfall_Timestamp(asynch,atoi(argv[3]),0);
	Asynch_Set_Last_Rainfall_Timestamp(asynch,atoi(argv[4]),0);

	//Set database to server 63 for rainfall data
	//sprintf(database_info,"dbname=ifloods host=s-iihr63.iihr.uiowa.edu port=5433 user=rainfall_feed password=r!Ain2012");
	//Asynch_Set_Database_Connection(asynch,database_info);

	//Advance solver for 1 hour = 60 min
	Asynch_Set_Total_Simulation_Time(asynch,60.0);
	Asynch_Advance(asynch,1);

	//Write current state to disk
	char timestamp[32];
	sprintf(timestamp,"%u",Asynch_Get_Last_Rainfall_Timestamp(asynch,0));
	Asynch_Take_System_Snapshot(asynch,timestamp);

	//Advance solver for another 23 hours
	Asynch_Set_Total_Simulation_Time(asynch,24.0*60.0);
	Asynch_Advance(asynch,1);

	//Create output files
	Asynch_Create_Output_Files(asynch);

	//Clean up
	Asynch_Free(asynch);
	return 0;
}

