#include <stdio.h>
#include <time.h>
#include <libpq-fe.h>
#include <string.h>
#include <unistd.h>
#include "asynch_interface.h"
#include "forecaster_methods.h"

int my_rank;
int np;

//!!!! The routines and structs here should probably go into a separate file !!!!

typedef struct CustomParams
{
	unsigned int ID;
	int offset;
} CustomParams;


int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user);
int Output_Timestamp(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user);
//double Output_Ratio(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user);
void Init_Output_User_forecastparams(asynchsolver* asynch);
void Free_Output_User_forecastparams(asynchsolver* asynch);
void Set_Output_User_forecastparams(asynchsolver* asynch,unsigned int offset);

void Init_Output_PeakflowUser_Offset(asynchsolver* asynch);
void Free_Output_PeakflowUser_Offset(asynchsolver* asynch);
void Set_Output_PeakflowUser_Offset(asynchsolver* asynch,unsigned int offset);

//Data assimilation functions
int AdjustDischarges_Scale(asynchsolver* asynch,ConnData* conninfo,unsigned int current_timestamp,double* areas,VEC** states);
int AdjustDischarges_Scale_ToSubsurface(asynchsolver* asynch,ConnData* conninfo,unsigned int current_timestamp,double* areas,VEC** states);
double* Find_Gauge_To_Use(asynchsolver* asynch);

int main(int argc,char* argv[])
{
	//Initialize MPI stuff
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	//Parse input
	if(argc < 3)
	{
		if(my_rank == 0)	printf("Command line parameter required:\nA universal variable file (.gbl),\nA forecast file (.fcst)\n");
		MPI_Finalize();
		return 1;
	}

	//Declare variables
	unsigned int i,j,k,current_offset;
	double total_time = 0.0;
	time_t start,start2,stop;
	asynchsolver* asynch;
	PGresult *res;
	MPI_Status status;
	char* query = (char*) malloc(1024*sizeof(char));
	Link* current;

	if(my_rank == 0)
		printf("\nBeginning initialization...\n*****************************\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Init asynch object and the river network
	asynch = Asynch_Init(MPI_COMM_WORLD);
	Asynch_Parse_GBL(asynch,argv[1]);
	Asynch_Load_System(asynch);

	//Load Forecast related data
	ForecastData* Forecaster = Init_ForecastData(argv[2],asynch->GlobalVars->string_size);
	if(!Forecaster)
		MPI_Abort(MPI_COMM_WORLD,1);

	//Setup output for link id, if needed
	int setup_id = Asynch_Check_Output(asynch,"LinkID");
	int setup_timestamp = Asynch_Check_Output(asynch,"Timestamp");
	//int setup_ratio = Asynch_Check_Output(asynch,"Ratio");
	//if( (setup_id || setup_timestamp || setup_ratio) != 0)
	if( (setup_id || setup_timestamp) != 0)
	{
		if(my_rank == 0)	printf("[%i]: Forecaster needs LinkID (%i), Timestamp (%i).\n",my_rank,setup_id,setup_timestamp);
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	Init_Output_User_forecastparams(asynch);
	Asynch_Set_Output(asynch,"LinkID",ASYNCH_INT,(void (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) &Output_Linkid,NULL,0);
	Asynch_Set_Output(asynch,"Timestamp",ASYNCH_INT,(void (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) &Output_Timestamp,NULL,0);
	//Asynch_Set_Output(asynch,"Ratio",ASYNCH_DOUBLE,(void (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) &Output_Ratio,NULL,0);

	Init_Output_PeakflowUser_Offset(asynch);

	//Get some values about the river system
	unsigned int N = Asynch_Get_Number_Links(asynch);
	unsigned int my_N = Asynch_Get_Local_Number_Links(asynch);
	char dump_filename[asynch->GlobalVars->string_size];

	//Create halt file
	CreateHaltFile(Forecaster->halt_filename);

	//Find the index of the forcing to use for forecasting
	unsigned int forecast_idx = Forecaster->forecasting_forcing;
	if(forecast_idx >= asynch->GlobalVars->num_forcings)
	{
		if(my_rank == 0)	printf("[%i]: Error: No forecasting forcing set.\n",my_rank);
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//Reserve space for backups
	VEC** backup = (VEC**) malloc(N*sizeof(VEC*));
	for(i=0;i<N;i++)
	{
		if(asynch->assignments[i] == my_rank || asynch->getting[i] == 1)
			backup[i] = v_get(asynch->GlobalVars->dim);
		else
			backup[i] = NULL;
	}

	if(my_rank == 0)
	{
		printf("\nModel type is %u.\nGlobal parameters are:\n",asynch->GlobalVars->type);
		Print_Vector(asynch->GlobalVars->global_params);
		printf("\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time = difftime(stop,start);
	if(my_rank == 0)	printf("Finished initializations. Total time: %f\n",total_time);
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(1);

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,my_N);
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Make the initial solve
	Asynch_Advance(asynch,0);

	//Stop the clock
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time += difftime(stop,start);

	//Output some data
	if(my_rank == 0)
	{
		printf("%i: The answer at ID %i at time %.12f is\n",my_rank,asynch->sys[asynch->my_sys[0]]->ID,asynch->sys[asynch->my_sys[0]]->last_t);
		Print_Vector(asynch->sys[asynch->my_sys[0]]->list->tail->y_approx);
		printf("Total time for calculations: %f\n",difftime(stop,start));
	}

	//Begin persistent calculations
	if(my_rank == 0)
		printf("\n\n===================================\nBeginning persistent calculations\n===================================\n");
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	//Make some initializations and checks
	//unsigned int history_time = 5*24*60*60;	//Amount of history to store for hydrographs and peakflows
	short unsigned int hr1 = 0;	//Hour of the day to perform maintainance on database
	short unsigned int hr2 = 12;	//Hour of the day to perform maintainance on database
	unsigned int wait_time = 120;	//Time to sleep if no rainfall data is available
	double forecast_time = 10.0*24*60;	//Time (mins) in future to make forecasts
	unsigned int num_tables = 10;
	//unsigned int num_rainsteps = 3;	//Number of rainfall intensities to use for the next forecast
	unsigned int num_rainsteps = Forecaster->num_rainsteps;	//Number of rainfall intensities to use for the next forecast
	//if(my_rank == 0 && asynch->GlobalVars->increment < num_rainsteps + 3)
	if(my_rank == 0 && asynch->forcings[forecast_idx]->increment < num_rainsteps + 3)
		printf("Warning: Increment for rain should probably be %u.\n",num_rainsteps + 3);
	asynch->forcings[forecast_idx]->increment = num_rainsteps;	//!!!! Not necessary, but makes me feel better. The solvers should really not do the last step where they download nothing. !!!!

	unsigned int nextraintime,repeat_for_errors,nextforcingtime;
	short int halt = 0;
	int isnull,error_code;
	int message_buffer[1 + asynch->GlobalVars->num_forcings];
	short int vac = 0;	//0 if no vacuum has occured, 1 if vacuum has occured (during a specific hour)
	unsigned int last_file = asynch->forcings[forecast_idx]->last_file;
	unsigned int first_file = asynch->forcings[forecast_idx]->first_file;
	k = 0;
	for(i=0;i<N;i++)
		if(backup[i] != NULL)	v_copy(asynch->sys[i]->list->tail->y_approx,backup[i]);

	double simulation_time_with_data = 0.0;
	//for(i=0;i<asynch->GlobalVars->num_forcings;i++)
	//	if(asynch->forcings[i]->flag == 5)	simulation_time_with_data = max(simulation_time_with_data,asynch->forcings[i]->file_time * asynch->forcings[i]->num_rainsteps);
	simulation_time_with_data = max(simulation_time_with_data,asynch->forcings[forecast_idx]->file_time * Forecaster->num_rainsteps);

	//Setup temp files
	Set_Output_User_forecastparams(asynch,first_file);
	Set_Output_PeakflowUser_Offset(asynch,first_file);
	Asynch_Set_Total_Simulation_Time(asynch,forecast_time);
	Asynch_Prepare_Temp_Files(asynch);

	//Prepare assimilation data
	double* areas = Find_Gauge_To_Use(asynch);
	ConnData* gauge_feed = CreateConnData("dbname=model_ifc host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");

	//Make some initializations to the database
	if(my_rank == 0)
	{
		printf("Making initializations to tables.\n");
		start = time(NULL);

		//Connect to hydrograph database
		ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

		//Make sure the hydrographs table exists
		sprintf(query,"SELECT 1 FROM pg_class WHERE relname='%s';",asynch->GlobalVars->hydro_table);
		res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
		if(!PQntuples(res))
		{
			PQclear(res);
			sprintf(query,"CREATE TABLE %s(link_id int,time int,ratio real,discharge real); ALTER TABLE %s SET (autovacuum_enabled = false, toast.autovacuum_enabled = false);",asynch->GlobalVars->hydro_table,asynch->GlobalVars->hydro_table);
			res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
			CheckResError(res,"creating hydrographs table");
		}
		else
		{
			PQclear(res);
			sprintf(query,"TRUNCATE %s;",asynch->GlobalVars->hydro_table);
			res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
			CheckResError(res,"truncating hydrographs table");
		}
		PQclear(res);

		//Make sure the hydrograph tables are set correctly
		CheckPartitionedTable(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->GlobalVars,Forecaster,num_tables,"archive_hydroforecast","forecast_time");

		//Clear the future hydrographs in archive
		DeleteFutureValues(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],num_tables,asynch->GlobalVars,"archive_hydroforecast",Forecaster->model_name,first_file,1);
/*
		sprintf(query,"DELETE FROM master_archive_hydroforecast_%s WHERE forecast_time >= %u;",Forecaster->model_name,first_file);
		res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
		CheckResError(res,"deleting future hydroforecasts");
		PQclear(res);
*/

		//Disconnect from hydrograph database, connect to peakflow database
		DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
		ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);

		//Make sure the peakforecast tables are set correctly
		//CheckPeakforecastTable(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT],asynch->GlobalVars,Forecaster,num_peakflow_tables);

		//Clear all future peakflows
		sprintf(query,"TRUNCATE %s;",asynch->GlobalVars->peak_table);
		res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]->conn,query);
		CheckResError(res,"truncating peakforecast table");
		PQclear(res);

		//Disconnect
		DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);

		stop = time(NULL);
		printf("Total time to initialize tables: %.2f.\n",difftime(stop,start));
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//Start the main loop
	while(!halt)
	{
		if(my_rank == 0)
		{
			time_t now = time(NULL);
			struct tm* now_info = localtime(&now);
			printf("\n\nPass %u\n",k);
			printf("Current time is %s",asctime(now_info));
		}

		//Clear buffers
		Flush_TransData(asynch->my_data);

		//Make some initializations
		//asynch->forcings[forecast_idx]->raindb_start_time = last_file;								//!!!! This all assumes one forcing from db !!!!
		first_file = last_file;
		last_file = last_file + (unsigned int) asynch->forcings[forecast_idx]->file_time * 60 * num_rainsteps;
		nextforcingtime = first_file + 60 * (unsigned int) rint(asynch->forcings[forecast_idx]->file_time) * (num_rainsteps-1);	//This is the actual timestamp of the last needed forcing data. This will be downloaded (unlike last_file)
		//nextforcingtime = first_file + 60 * (unsigned int) rint(asynch->forcings[forecast_idx]->file_time) * num_rainsteps;

		//Assimilate the observations
		//AdjustDischarges_Scale(asynch,gauge_feed,first_file,areas,backup);
		AdjustDischarges_Scale_ToSubsurface(asynch,gauge_feed,first_file,areas,backup);

		//Reset each link
		Asynch_Set_System_State(asynch,0.0,backup);
		Set_Output_User_forecastparams(asynch,first_file);
		Set_Output_PeakflowUser_Offset(asynch,first_file);
		Asynch_Write_Current_Step(asynch);
		Asynch_Set_Forcing_State(asynch,forecast_idx,0.0,first_file,last_file);

		for(i=0;i<asynch->GlobalVars->num_forcings;i++)	//Set any other database forcings to begin at first_file
		{
			if(asynch->forcings[i]->flag == 3)
				Asynch_Set_Forcing_State(asynch,i,0.0,first_file,asynch->forcings[i]->last_file);
		}

		//Check if a vacuum should be done
		//This will happen at hr1
		if(my_rank == 0)	PerformTableMaintainance(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->GlobalVars,Forecaster,&vac,hr1,num_tables,"archive_hydroforecast");

		//Make sure all buffer flushing is done
		MPI_Barrier(MPI_COMM_WORLD);

		//Dump data for debugging and recovery
		if(k % 96 == 0)
		{
			//sprintf(dump_filename,"%u",last_file);
			sprintf(dump_filename,"/%u",first_file);
			Asynch_Take_System_Snapshot(asynch,dump_filename);
		}

		//Find the next time where rainfall occurs
		do
		{
			if(my_rank == 0)
			{
				//Init message_buffer
				//message_buffer[0] = 1;

				//Check if connection to SQL database is still good
				//for(i=0;i<asynch->GlobalVars->num_forcings;i++)
				{
					//message_buffer[i+1] = 0x7FFFFFFF;	//Largest 32 bit int	!!!! I don't think this is necessary !!!!

					//if(asynch->forcings[i]->flag == 5)
					{
						ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_FORCING_START + forecast_idx]);

						//Find the next rainfall time
						time(&start);
						sprintf(query,"SELECT min(unix_time) FROM rain_maps5_index WHERE unix_time >= %u AND link_count > -1;",nextforcingtime);	//!!!! Should vary with forcing (i) !!!!
						res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_FORCING_START + forecast_idx]->conn,query);
						CheckResError(res,"checking for new rainfall data");
						//time(&stop);
						//printf("Total time to check for new rainfall data: %f.\n",difftime(stop,start));

						isnull = PQgetisnull(res,0,0);
/*
						if(!isnull)
						{
							message_buffer[0] = 0;
							int value = atoi(PQgetvalue(res,0,0));
							message_buffer[i+1] = (value < message_buffer[i+1]) ? value : message_buffer[i+1];
							PQclear(res);
							DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_FORCING_START + i]);
						}
						else
						{
							message_buffer[0] = 1;
							message_buffer[i+1] = 0;
							PQclear(res);
							DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_FORCING_START + i]);
							break;
						}
*/
						PQclear(res);
						DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_FORCING_START + forecast_idx]);

						//Check if there is new gauge data
						ConnectPGDB(gauge_feed);

						sprintf(query,"SELECT min(unix_time) FROM _data_usgs WHERE unix_time >= %u;",nextforcingtime);
						res = PQexec(gauge_feed->conn,query);
						CheckResError(res,"checking for new gauge data");
						isnull |= PQgetisnull(res,0,0);
						PQclear(res);

						//Check if there is enough gauge data
						unsigned int num_gauges;
						sprintf(query,"SELECT count(DISTINCT ifis_id) FROM _data_usgs WHERE unix_time >= %u;",nextforcingtime);
						res = PQexec(gauge_feed->conn,query);
						error_code = CheckResError(res,"counting number of gauges");
						num_gauges = (error_code) ? 0 : atoi(PQgetvalue(res,0,0));

printf("[%i]: Got data for %u gauges.\n",my_rank,num_gauges);

						if(num_gauges < 140)	isnull = 1;
						PQclear(res);

						DisconnectPGDB(gauge_feed);

						time(&stop);
						printf("Total time to check for new rainfall and gauge data: %f.\n",difftime(stop,start));
					}
				}
			}
			MPI_Bcast(&isnull,1,MPI_INT,0,MPI_COMM_WORLD);
			//MPI_Bcast(message_buffer,1+asynch->GlobalVars->num_forcings,MPI_INT,0,MPI_COMM_WORLD);
			//isnull = message_buffer[0];

			//if(!isnull)	//!!!! Seems I don't need nextraintime. Do I need to even send more than just isnull? !!!!
				//nextraintime = message_buffer[1];
			//else
			if(isnull)
			{
				if(my_rank == 0)
				{
					printf("No value returned from SQL database for forcing or gauges. %u %u\n",last_file,isnull);
					//PerformTableMaintainance(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->GlobalVars,Forecaster,&vac,hr1,hr2,first_file - history_time,num_tables);
					PerformTableMaintainance(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->GlobalVars,Forecaster,&vac,hr1,num_tables,"archive_hydroforecast");
				}

				halt = CheckFinished(Forecaster->halt_filename);
				if(halt)
				{
					sprintf(dump_filename,"/%u",first_file);
					Asynch_Take_System_Snapshot(asynch,dump_filename);
					//DataDump(sys,N,assignments,GlobalVars,first_file);
				}
				else
				{
					fflush(stdout);
					sleep(wait_time);
				}
			}
		} while(isnull && !halt);

		if(halt)	break;

		//Read in next set of rainfall data

		//Initialize some data for the first phase of calculations
		//GlobalVars->maxtime = GlobalVars->file_time * num_rainsteps;
		Asynch_Set_Total_Simulation_Time(asynch,simulation_time_with_data);		// !!!! This may not work for multiple forcings for forecasting. How do you handle different time resolutions? !!!!
		//conninfo->time_offset = first_file;
		current_offset = first_file;
		Set_Output_User_forecastparams(asynch,current_offset);
		Set_Output_PeakflowUser_Offset(asynch,current_offset);

		MPI_Barrier(MPI_COMM_WORLD);
		time(&start);
if(my_rank == 0)
printf("first: %u last: %u\n",first_file,last_file);
		Asynch_Advance(asynch,1);

		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			time(&stop);
			printf("Time for first phase calculations: %.2f\n",difftime(stop,start));
		}

		//Flush communication buffers
		Flush_TransData(asynch->my_data);

		//Reset the links (mostly) and make a backup for the second phase
		for(i=0;i<N;i++)	//Set time to 0.0
		{
			current = asynch->sys[i];
			if(current->list != NULL)
			{
				while(current->current_iterations > 1)
				{
					Remove_Head_Node(current->list);
					(current->current_iterations)--;
				}

				current->steps_on_diff_proc = 1;
				current->iters_removed = 0;
				current->rejected = 0;
				if(current->numparents == 0)	current->ready = 1;
				else				current->ready = 0;
//				v_copy(current->list->head->y_approx,backup[i]);
			}
		}

		//Make second phase calculations
		MPI_Barrier(MPI_COMM_WORLD);
		time(&start);
		//GlobalVars->maxtime = forecast_time;
		Asynch_Set_Total_Simulation_Time(asynch,forecast_time);
		Asynch_Deactivate_Forcing(asynch,forecast_idx);
		Asynch_Advance(asynch,1);
		Asynch_Activate_Forcing(asynch,forecast_idx);
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			time(&stop);
			printf("Time for second phase calculations: %.2f\n",difftime(stop,start));
		}

		//Output some data
		if(my_rank == 0)
		{
			printf("[%i]: The answer at ID %i at time %.12f is\n",my_rank,asynch->sys[asynch->my_sys[0]]->ID,asynch->sys[asynch->my_sys[0]]->last_t);
			Print_Vector(asynch->sys[asynch->my_sys[0]]->list->tail->y_approx);
		}

		//Upload the peak data to the database **********************************************************************************************

		//Adjust the table peak
		MPI_Barrier(MPI_COMM_WORLD);
		start = time(NULL);

		repeat_for_errors = Asynch_Create_Peakflows_Output(asynch);
		while(repeat_for_errors > 0)
		{
			if(my_rank == 0)	printf("[%i]: Attempting resend of peakflow data.\n",my_rank);
			sleep(5);
			repeat_for_errors = Asynch_Create_Peakflows_Output(asynch);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			stop = time(NULL);
			printf("[%i]: Total time to transfer peak flow data: %.2f\n",my_rank,difftime(stop,start));
		}

		//Upload the hydrographs to the database ********************************************************************************************
		MPI_Barrier(MPI_COMM_WORLD);
		start = time(NULL);

		//Adjust the table hydrographs
		if(my_rank == 0)
		{
			//Make sure database connection is still good
			ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

			sprintf(query,"TRUNCATE %s;",asynch->GlobalVars->hydro_table);
			res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
			CheckResError(res,"deleting hydrographs");
			PQclear(res);

			DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		repeat_for_errors = Asynch_Create_Output(asynch,NULL);
		while(repeat_for_errors > 0)
		{
			if(my_rank == 0)	printf("[%i]: Attempting resend of hydrographs data.\n",my_rank);
			sleep(5);
			repeat_for_errors = Asynch_Create_Output(asynch,NULL);
		}

		//Call functions
		if(my_rank == 0)
		{
			//Connect to database
			ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

			//Functions for displaying data on IFIS
			if(Forecaster->ifis_display)
			{
				//Stage
				repeat_for_errors = 1;
				while(repeat_for_errors)
				{
					repeat_for_errors = 0;
					sprintf(query,"SELECT get_stages_%s();",Forecaster->model_name);
					res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
					repeat_for_errors = repeat_for_errors || CheckResError(res,"calling stage function");
					PQclear(res);
					if(repeat_for_errors)
					{
						printf("[%i]: Attempting to call stage function again...\n",my_rank);
						sleep(5);
						CheckConnConnection(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
					}
				}

				//Warnings
				repeat_for_errors = 1;
				while(repeat_for_errors)
				{
					repeat_for_errors = 0;
					sprintf(query,"SELECT update_warnings_%s();",Forecaster->model_name);
					res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
					repeat_for_errors = repeat_for_errors || CheckResError(res,"calling warnings function");
					PQclear(res);
					if(repeat_for_errors)
					{
						printf("[%i]: Attempting to call warning function again...\n",my_rank);
						sleep(5);
						CheckConnConnection(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
					}
				}
			}

			//Stage archive
			repeat_for_errors = 1;
			while(repeat_for_errors)
			{
				repeat_for_errors = 0;
				sprintf(query,"ALTER TABLE master_archive_hydroforecast_%s ALTER COLUMN forecast_time SET DEFAULT %u;",Forecaster->model_name,current_offset);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"setting default value");
				PQclear(res);

				sprintf(query,"SELECT copy_to_archive_hydroforecast_%s();",Forecaster->model_name);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"calling stage archive function");
				PQclear(res);

				sprintf(query,"ALTER TABLE master_archive_hydroforecast_%s ALTER COLUMN forecast_time DROP DEFAULT;",Forecaster->model_name);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"dropping default value");
				PQclear(res);

				if(repeat_for_errors)
				{
					printf("[%i]: Attempting to call stage archive function again...\n",my_rank);
					sleep(5);
					CheckConnConnection(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
				}
			}

			//Disconnect
			DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
		}

		if(my_rank == 0)
		{
			time(&stop);
			printf("[%i]: Total time to transfer hydrograph data: %.2f\n",my_rank,difftime(stop,start));
		}
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);

		//Check if program has received a terminate signal **********************************************************************************
		k++;
		halt = CheckFinished(Forecaster->halt_filename);

		//If stopping, make a .rec file
		if(halt)
		{
			for(i=0;i<N;i++)
			{
				current = asynch->sys[i];
				if(current->list != NULL)
					v_copy(backup[i],current->list->tail->y_approx);
			}

			sprintf(dump_filename,"/%u",last_file);
			Asynch_Take_System_Snapshot(asynch,dump_filename);
		}
	}

	//Clean up **********************************************************************************************************************************
	free(query);
	for(i=0;i<N;i++)	v_free(backup[i]);
	free(backup);
	Free_ForecastData(&Forecaster);
	Asynch_Delete_Temporary_Files(asynch);
	Free_Output_PeakflowUser_Offset(asynch);
	Free_Output_User_forecastparams(asynch);
	Asynch_Free(asynch);
	free(areas);
	ConnData_Free(gauge_feed);
	return 0;
}


//Output functions ****************************************************************************
int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user)
{
	CustomParams* forecastparams = (CustomParams*) user;
	return forecastparams->ID;
}

int Output_Timestamp(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user)
{
	CustomParams* forecastparams = (CustomParams*) user;
	return (int)(round(t * 60.0 + forecastparams->offset) + 0.1);
}


//Custom parameters for forecasting ***********************************************************
void Init_Output_User_forecastparams(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		sys[my_sys[i]]->output_user = malloc(sizeof(CustomParams));
}

void Free_Output_User_forecastparams(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
	{
		free(sys[my_sys[i]]->output_user);
		sys[my_sys[i]]->output_user = NULL;
	}
}

void Set_Output_User_forecastparams(asynchsolver* asynch,unsigned int offset)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;
	CustomParams* forecastparams;

	for(i=0;i<my_N;i++)
	{
		forecastparams = (CustomParams*) sys[my_sys[i]]->output_user;
		forecastparams->ID = sys[my_sys[i]]->ID;
		forecastparams->offset = offset;
		//forecastparams->Q_TM = sys[my_sys[i]]->list->tail->y_approx->ve[0];
	}
}

void Init_Output_PeakflowUser_Offset(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		sys[my_sys[i]]->peakoutput_user = malloc(sizeof(unsigned int));
}

void Free_Output_PeakflowUser_Offset(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
	{
		free(sys[my_sys[i]]->peakoutput_user);
		sys[my_sys[i]]->peakoutput_user = NULL;
	}
}

void Set_Output_PeakflowUser_Offset(asynchsolver* asynch,unsigned int offset)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		*(unsigned int*)(sys[my_sys[i]]->peakoutput_user) = offset;
}




//Downloads gauge data and scales discharge upstream by area
//This assumes the states array has been set for all non-discharge states (or will be set)
//!!!! This assumes the number gauges could vary. If not needed, calculate gauge dependence in an init routine. !!!!
int AdjustDischarges_Scale(asynchsolver* asynch,ConnData* conninfo,unsigned int current_timestamp,double* areas,VEC** states)
{
	//Unpack
	UnivVars* GlobalVars = asynch->GlobalVars;
	unsigned int N = asynch->N,*my_sys = asynch->my_sys,my_N = asynch->my_N,**id_to_loc = asynch->id_to_loc;
	Link** sys = asynch->sys;
	short int *getting = asynch->getting;

	unsigned int i,j,num_gauges,loc,*upstream = (unsigned int*) malloc(N*sizeof(unsigned int)),num_upstream,prev_loc,curr_loc,*ids;
	short int *locs_set = (short int*) calloc(N,sizeof(short int));
	double *locs_newq = (double*) malloc(N*sizeof(double)),*qs;
	unsigned int *counter = (unsigned int*) malloc(N*sizeof(unsigned int));
	PGresult *res;
	char query[GlobalVars->query_size];
	Link* current;
	double ratio;

	for(i=0;i<N;i++)
	{
		locs_newq[i] = -1.0;
	}

	//Get data from database
	if(my_rank == 0)
	{
		ConnectPGDB(conninfo);
		sprintf(query,"SELECT B.link_id,A.discharge*0.0283168 FROM _data_usgs AS A,_pois_adv AS B WHERE A.ifis_id = B.id AND unix_time = %u ORDER BY B.link_id;",current_timestamp);
printf("Query: %s\n",query);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"downloading gauge readings");
		num_gauges = PQntuples(res);
printf("Got %u gauges.\n",num_gauges);

		//Separate data from res
		MPI_Bcast(&num_gauges,1,MPI_INT,0,MPI_COMM_WORLD);
		ids = (unsigned int*) malloc(num_gauges*sizeof(unsigned int));
		qs = (double*) malloc(num_gauges*sizeof(double));
		for(i=0;i<num_gauges;i++)
		{
			ids[i] = atoi(PQgetvalue(res,i,0));
			qs[i] = atof(PQgetvalue(res,i,1));
		}

		PQclear(res);
		DisconnectPGDB(conninfo);
	}
	else
	{
		MPI_Bcast(&num_gauges,1,MPI_INT,0,MPI_COMM_WORLD);
		ids = (unsigned int*) malloc(num_gauges*sizeof(unsigned int));
		qs = (double*) malloc(num_gauges*sizeof(double));
	}
	MPI_Bcast(ids,num_gauges,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(qs,num_gauges,MPI_DOUBLE,0,MPI_COMM_WORLD);

	//Store gauge readings
	for(i=0;i<num_gauges;i++)
	{
		loc = find_link_by_idtoloc(ids[i],id_to_loc,N);
		if(loc < N)
		{
			locs_set[loc] = 1;
			locs_newq[loc] = qs[i];
		}
	}

	//Set the discharges downstream from each gauge. This is for locations with no downstream gauages.
	//This follows each gauge downstream until a link is found with locs_set or an outlet.
	//counter is set to help track the closest gauge.
	for(i=0;i<N;i++)	counter[i] = N;

	for(i=0;i<num_gauges;i++)
	{
		loc = find_link_by_idtoloc(ids[i],id_to_loc,N);
		if(loc < N)
		{
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
						locs_newq[curr_loc] = areas[current->location] * locs_newq[loc] / areas[loc];
					}
					prev_loc = curr_loc;
					current = current->c;
					if(current)	curr_loc = current->location;
					else		break;
				}
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
				ratio = locs_newq[current->location] / areas[current->location];
				for(j=0;j<num_upstream;j++)
				{
					locs_set[upstream[j]] = 1;
					locs_newq[upstream[j]] = areas[upstream[j]] * ratio;
				}
			}
		}
	}

	//Set the determined discharge. If a link's discharge was not determined by the above process, then it lies in a totally ungauged basin.
	for(i=0;i<N;i++)
	{
		if(states[i] && locs_set[i])
			states[i]->ve[0] = locs_newq[i];
	}

	//Clean up
	free(locs_set);
	free(upstream);
	free(locs_newq);
	free(counter);
	free(qs);
	free(ids);
	return 0;
}


//Downloads gauge data and scales discharge upstream by area. Water is put or pulled from the subsurface to preserve volume.
//This assumes the states array has been set for all states.
//This requires model 252 or 253
//!!!! This assumes the number gauges could vary. If not needed, calculate gauge dependence in an init routine. !!!!
int AdjustDischarges_Scale_ToSubsurface(asynchsolver* asynch,ConnData* conninfo,unsigned int current_timestamp,double* areas,VEC** states)
{
	//Unpack
	UnivVars* GlobalVars = asynch->GlobalVars;
	unsigned int N = asynch->N,*my_sys = asynch->my_sys,my_N = asynch->my_N,**id_to_loc = asynch->id_to_loc;
	Link** sys = asynch->sys;
	short int *getting = asynch->getting;

	unsigned int i,j,num_gauges,loc,*upstream = (unsigned int*) malloc(N*sizeof(unsigned int)),num_upstream,prev_loc,curr_loc,*ids;
	short int *locs_set = (short int*) calloc(N,sizeof(short int));
	double *locs_newq = (double*) malloc(N*sizeof(double)),*qs;
	unsigned int *counter = (unsigned int*) malloc(N*sizeof(unsigned int));
	PGresult *res;
	char query[GlobalVars->query_size];
	Link* current;
	double ratio;

	for(i=0;i<N;i++)
	{
		locs_newq[i] = -1.0;
	}

	//Get data from database
	if(my_rank == 0)
	{
		ConnectPGDB(conninfo);
		sprintf(query,"SELECT B.link_id,A.discharge*0.0283168 FROM _data_usgs AS A,_pois_adv AS B WHERE A.ifis_id = B.id AND unix_time = %u ORDER BY B.link_id;",current_timestamp);
printf("Query: %s\n",query);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"downloading gauge readings");
		num_gauges = PQntuples(res);
printf("Got %u gauges.\n",num_gauges);

		//Separate data from res
		MPI_Bcast(&num_gauges,1,MPI_INT,0,MPI_COMM_WORLD);
		ids = (unsigned int*) malloc(num_gauges*sizeof(unsigned int));
		qs = (double*) malloc(num_gauges*sizeof(double));
		for(i=0;i<num_gauges;i++)
		{
			ids[i] = atoi(PQgetvalue(res,i,0));
			qs[i] = atof(PQgetvalue(res,i,1));
		}

		PQclear(res);
		DisconnectPGDB(conninfo);
	}
	else
	{
		MPI_Bcast(&num_gauges,1,MPI_INT,0,MPI_COMM_WORLD);
		ids = (unsigned int*) malloc(num_gauges*sizeof(unsigned int));
		qs = (double*) malloc(num_gauges*sizeof(double));
	}
	MPI_Bcast(ids,num_gauges,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(qs,num_gauges,MPI_DOUBLE,0,MPI_COMM_WORLD);

	//Store gauge readings
	for(i=0;i<num_gauges;i++)
	{
		loc = find_link_by_idtoloc(ids[i],id_to_loc,N);
		if(loc < N)
		{
			locs_set[loc] = 1;
			locs_newq[loc] = qs[i];
		}
	}

	//Set the discharges downstream from each gauge. This is for locations with no downstream gauages.
	//This follows each gauge downstream until a link is found with locs_set or an outlet.
	//counter is set to help track the closest gauge.
	for(i=0;i<N;i++)	counter[i] = N;

	for(i=0;i<num_gauges;i++)
	{
		loc = find_link_by_idtoloc(ids[i],id_to_loc,N);
		if(loc < N)
		{
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
						locs_newq[curr_loc] = areas[current->location] * locs_newq[loc] / areas[loc];
					}
					prev_loc = curr_loc;
					current = current->c;
					if(current)	curr_loc = current->location;
					else		break;
				}
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
				ratio = locs_newq[current->location] / areas[current->location];
				for(j=0;j<num_upstream;j++)
				{
					locs_set[upstream[j]] = 1;
					locs_newq[upstream[j]] = areas[upstream[j]] * ratio;
				}
			}
		}
	}

	//Set the determined discharge. If a link's discharge was not determined by the above process, then it lies in a totally ungauged basin.
	double one_m_lambda_1 = 1.0 - asynch->GlobalVars->global_params->ve[1],tau,S_start,S_final,S_diff,S_avail,A_h;
	for(i=0;i<N;i++)
	{
		if(states[i] && locs_set[i])
		{
			tau = 1.0 / sys[i]->params->ve[3];
			A_h = sys[i]->params->ve[2];

			//Calculate how much water to pull from the subsurface
			S_start = pow(states[i]->ve[0] * tau / 60.0,one_m_lambda_1);
			S_final = pow(locs_newq[i] * tau / 60.0,one_m_lambda_1);
			S_diff = S_final - S_start;

			//Calculate how much water is available
			S_avail = A_h * states[i]->ve[3];

			if(S_diff < 0.0)	//Add water to hillslope
			{
				states[i]->ve[0] = locs_newq[i];
				states[i]->ve[3] += S_diff / A_h;
			}
			else		//Remove water from hillslope
			{
				if(S_avail >= S_diff)	//Got enough water
				{
					states[i]->ve[0] = locs_newq[i];
					states[i]->ve[3] += S_diff / A_h;
				}
				else	//Not enough
				{
					states[i]->ve[0] = locs_newq[i];
					states[i]->ve[3] = 0.0;
				}
			}
		}
	}

	//Clean up
	free(locs_set);
	free(upstream);
	free(locs_newq);
	free(counter);
	free(qs);
	free(ids);
	return 0;
}



//Provides a list of areas for all links to all procs
double* Find_Gauge_To_Use(asynchsolver* asynch)
{
	//Unpack
	unsigned int N = asynch->N,*assignments = asynch->assignments,area_idx = asynch->GlobalVars->area_idx;
	Link** sys = asynch->sys;

	unsigned int i;
	double* areas = (double*) malloc(N*sizeof(double));

	//Broadcast areas
	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank)
		{
			MPI_Bcast(&(sys[i]->params->ve[area_idx]),1,MPI_DOUBLE,my_rank,MPI_COMM_WORLD);
			areas[i] = sys[i]->params->ve[area_idx];
		}
		else
			MPI_Bcast(&(areas[i]),1,MPI_DOUBLE,assignments[i],MPI_COMM_WORLD);
	}

	return areas;
}


