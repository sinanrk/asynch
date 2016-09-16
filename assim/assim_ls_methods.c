#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#if !defined(_MSC_VER)
#define ASYNCH_SLEEP sleep
#else
#include <windows.h>
#define ASYNCH_SLEEP Sleep
#endif

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#include "assim_ls_methods.h"

//!!!! Use interface instead !!!!
void ResetSysLS(Link** sys,unsigned int N,UnivVars* GlobalVars,double t_0,double* backup,unsigned int problem_dim,unsigned int num_forcings,TransData* my_data)
{
	unsigned i,j,k,l;
	Link* current;

	Flush_TransData(my_data);

	for(i=0;i<N;i++)
	{
		current = sys[i];
		if(current->list != NULL)
		{
			while(current->current_iterations > 1)
			{
				Remove_Head_Node(current->list);
				(current->current_iterations)--;
			}
			current->list->head->t = t_0;
			current->last_t = t_0;
			current->steps_on_diff_proc = 1;
			current->iters_removed = 0;
			current->rejected = 0;
			if(current->numparents == 0)	current->ready = 1;
			else				current->ready = 0;
			for(j=0;j<problem_dim;j++)	current->list->head->y_approx.ve[j] = backup[i*problem_dim+j];
			//v_copy(backup[i],current->list->head->y_approx);

			//Reset the next_save time
			if(current->save_flag)
			{
				current->next_save = t_0;		//!!!! This forces the print times to match up with the assimilation times !!!!
				//current->disk_iterations = 1;
			}

			//Reset peak flow information
			current->peak_time = t_0;
			v_copy(current->list->head->y_approx,current->peak_value);

			//Set hydrograph scale
			//current->Q_TM = backup[i]->ve[0];

			//Reset current state
			if(current->state_check != NULL)
				current->state = current->state_check(current->list->head->y_approx,GlobalVars->global_params,current->params,current->qvs,current->dam);
			current->list->head->state = current->state;

			//Set rainfall
			if(current->forcing_buff)
			{
				for(k=0;k<num_forcings;k++)
				{
					if(current->forcing_buff[k])
					{
						//Find the right index in rainfall
						for(l=0;l<current->forcing_buff[k]->n_times-1;l++)
							if( current->forcing_buff[k]->rainfall[l][0] <= t_0 && t_0 < current->forcing_buff[k]->rainfall[l+1][0] )	break;
						double rainfall_buffer = current->forcing_buff[k]->rainfall[l][1];
						current->forcing_values[k] = rainfall_buffer;
						current->forcing_indices[k] = l;

						//Find and set the new change in rainfall
						for(j=l+1;j<current->forcing_buff[k]->n_times;j++)
						{
							if( fabs(current->forcing_buff[k]->rainfall[j][1] - rainfall_buffer) > 1e-12 )
							{
								current->forcing_change_times[k] = current->forcing_buff[k]->rainfall[j][0];
								break;
							}
						}
						if(j == current->forcing_buff[k]->n_times)
							current->forcing_change_times[k] = current->forcing_buff[k]->rainfall[j-1][0];

						//Select new step size
						//current->h = InitialStepSize(current->last_t,current,GlobalVars,workspace);
					}
				}
			}
		}
	}
}

//Finds the link ids upstream from every link in data_locs. If trim is 1, then only links which can affect the links in data_locs (assuming a constant channel velocity) are used.
void Find_Upstream_Links(asynchsolver* asynch,unsigned int problem_dim,short int trim,double inc,unsigned int steps_to_use,unsigned int* data_locs,unsigned int numdata)
{
	Link **sys = asynch->sys,*current;
	unsigned int N = asynch->N,parentsval,leaves_size = 0,i,j,l,m,**id_to_loc = asynch->id_to_loc;
	int *assignments = asynch->assignments;
	UnivVars *GlobalVars = asynch->GlobalVars;
	Link **leaves = (Link**) malloc(N*sizeof(Link*));
	Link **stack = (Link**) malloc(N*sizeof(Link*));
	unsigned short int* getting = asynch->getting;
	upstream_data* updata;

	//Find leaves
	for(i=0;i<N;i++)
		if(sys[i]->numparents == 0)	leaves[leaves_size++] = sys[i];

	unsigned int* temp_numupstream = (unsigned int*) calloc(N,sizeof(unsigned int));
	for(i=0;i<leaves_size;i++)
		temp_numupstream[leaves[i]->location] = 1;

	//Count upstream links
	for(i=0;i<leaves_size;i++)
	{
		for(current = leaves[i]->c; current; current = current->c)
		{
			parentsval = 0;
			for(j=0;j<current->numparents;j++)	parentsval += (temp_numupstream[current->parents[j]->location] > 0);

			if(parentsval == current->numparents)	//All parents have temp_numupstream set
			{
				temp_numupstream[current->location] = 1;
				for(j=0;j<current->numparents;j++)
					temp_numupstream[current->location] += temp_numupstream[current->parents[j]->location];
			}
			else
				break;
		}
	}

	//Set the upstream links
	unsigned int** temp_upstream = (unsigned int**) malloc(N*sizeof(unsigned int*));	//temp_upstream[i] is a list of all links upstream from link i
	for(i=0;i<N;i++)
		temp_upstream[i] = (unsigned int*) malloc(temp_numupstream[sys[i]->location] * sizeof(unsigned int));
	unsigned int* counter = (unsigned int*) calloc(N,sizeof(unsigned int));

	unsigned int stack_size = leaves_size;
	for(i=0;i<leaves_size;i++)	stack[i] = leaves[i];

	while(stack_size > 0)
	{
		current = stack[stack_size-1];
		l = current->location;

		//Add this link to its own upstream list
		temp_upstream[l][counter[l]] = l;
		counter[l]++;

		//Add each parents' upstream list
		for(i=0;i<current->numparents;i++)
		{
			m = current->parents[i]->location;
			if(sys[m]->res)	continue;
			for(j=0;j<counter[m];j++)
				temp_upstream[l][counter[l]+j] = temp_upstream[m][j];
			counter[l] += counter[m];
		}

		stack_size--;

		//If every parent of current's child has an upstream list determined, add it to the stack
		if(current->c != NULL)
		{
			parentsval = 0;
			for(i=0;i<current->c->numparents;i++)
			{
				m = current->c->parents[i]->location;
				parentsval += (counter[m] > 0);
			}

			if(parentsval == current->c->numparents)
			{
				stack[stack_size] = current->c;
				stack_size++;
			}
		}
	}

	//Move the data from temp_upstream into the child upstream
	short int* used = (short int*) calloc(N,sizeof(short int));	//1 if temp_upstream[i] was used, 0 if not
	for(i=0;i<N;i++)
	{
		//if(assignments[i] == my_rank || getting[i])
		{
			sys[i]->user = malloc(sizeof(upstream_data));
			updata = (upstream_data*) (sys[i]->user);

			updata->upstream = (unsigned int**) malloc(sys[i]->numparents * sizeof(unsigned int*));
			updata->num_upstream = (unsigned int*) malloc(sys[i]->numparents * sizeof(unsigned int));
			updata->fit_states = NULL;
			updata->fit_to_universal = NULL;
			for(j=0;j<sys[i]->numparents;j++)
			{
				updata->upstream[j] = temp_upstream[sys[i]->parents[j]->location];
				updata->num_upstream[j] = temp_numupstream[sys[i]->parents[j]->location];
				used[sys[i]->parents[j]->location] = 1;
			}
/*
			//Calculate the dimension at each link
			updata->dim = problem_dim + problem_dim + (problem_dim-1)*(problem_dim-1);	//Model eqs + variational eqs from this link
			for(j=0;j<sys[i]->numparents;j++)
				updata->dim += updata->num_upstream[j] * problem_dim;	//Variational eqs from upstream !!!! Too high? !!!!
*/
		}
	}


	//Cleanup
	for(i=0;i<N;i++)
	{
		if(!used[i])	free(temp_upstream[i]);
	}
	free(temp_upstream);
	free(temp_numupstream);
	free(counter);
	free(leaves);
	free(stack);
	free(used);

	//Remove extra links from the upstream lists
	//!!!! This needs to be generalized for the case of many outlets !!!!
	if(trim)	//!!!! Does this work correctly for the entire state? !!!!
	{
		char query[1028];
		unsigned int loc;
		PGresult *res;
		double *distance = (double*) malloc(N*sizeof(double)),influence_radius;
		double speed = 3.0 * 60.0;	//In m/min

		//!!!! Hard coding right now. Blah... !!!!
		unsigned int outlet = asynch->GlobalVars->outletlink;
		//unsigned int outlet = 434478;  //Turkey River above French Hollow
		//unsigned int outlet = 434514;	//Turkey River at Garber
		//unsigned int outlet = 307864;  //Half Squaw Creek
		//unsigned int outlet = 292254;	//Squaw Creek at Ames

		if(my_rank == 0)
		{
			//ConnData* conninfo = CreateConnData("dbname=arch_usgs host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
			ConnData* conninfo = CreateConnData("dbname=model_assim host=s-iihr51.iihr.uiowa.edu port=5432 user=automated_solver password=C5.pfest0");
			ConnectPGDB(conninfo);

			sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id FROM env_master_km AS nodeX, env_master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) \
					SELECT distinct env_master_km.link_id,to_border*1609.34 FROM env_master_km,subbasin WHERE subbasin.link_id = env_master_km.link_id;",outlet);
			//sprintf(query,"SELECT distinct env_master_km.link_id,to_border*1609.34 FROM env_master_km;");
			res = PQexec(conninfo->conn,query);
			if(CheckResError(res,"getting list of distances to outlet"))
				MPI_Abort(MPI_COMM_WORLD,1);
			i = PQntuples(res);
			if(i != N)
				printf("Error: got a different number of links for the distances to outlet than links in network. (%u vs %u)\n",i,N);
			else
			{
				//Sort the data
				for(i=0;i<N;i++)
				{
					loc = find_link_by_idtoloc(atoi(PQgetvalue(res,i,0)),id_to_loc,N);
					if(loc >= N)
					{
						i = loc;
						break;
					}
					distance[loc] = atof(PQgetvalue(res,i,1));
				}
			}

			//Clean up db connection
			PQclear(res);
			DisconnectPGDB(conninfo);
			ConnData_Free(conninfo);

//printf("!!!! Loading distances for test basin !!!!\n");
//distance[0] = 0;distance[1] = 1;distance[2] = 2;distance[3] = 2;distance[4] = 3;distance[5] = 1;distance[6] = 3;distance[7] = 4;distance[8] = 2;distance[9] = 3;distance[10] = 3;
		}
//speed = 0.021;
		MPI_Bcast(&i,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		if(i != N)	return;
		MPI_Bcast(distance,N,MPI_DOUBLE,0,MPI_COMM_WORLD);

		//Calculate the radius of influence
		influence_radius = speed * steps_to_use * inc;

		//For each link, find the gauges influenced. Then, check which upstream links influence one of those gauges.
		unsigned int *influenced_gauges = (unsigned int*) malloc(numdata*sizeof(unsigned int)),num_influenced,drop;
		double difference;
		for(i=0;i<N;i++)
		{
			current = sys[i];
			updata = (upstream_data*) current->user;

			//Find influenced gauges
			num_influenced = 0;
			for(j=0;j<numdata;j++)
			{
				difference = distance[i] - distance[data_locs[j]];
				if(-1e-12 < difference && difference < influence_radius)	//!!!! Does this work if the network is disconnected? !!!!
{
					influenced_gauges[num_influenced++] = data_locs[j];
//if(sys[data_locs[j]]->ID == 399711)
//printf("%u\n",current->ID);
}
			}

			//Find which upstream links are also influenced
			for(j=0;j<current->numparents;j++)
			{
				drop = 0;
				for(l=0;l<updata->num_upstream[j];l++)
				{
					for(m=0;m<num_influenced;m++)
					{
						if( distance[updata->upstream[j][l]] - distance[influenced_gauges[m]] < influence_radius )
							break;
					}

					if(m == num_influenced)	//Upstream link does not influence any gauged location
						drop++;
					else
						updata->upstream[j][l-drop] = updata->upstream[j][l];
				}
				updata->num_upstream[j] -= drop;
				updata->upstream[j] = (unsigned int*) realloc(updata->upstream[j],updata->num_upstream[j]*sizeof(unsigned int));
			}
		}


/*
for(i=0;i<N;i++)
{
current = sys[i];
updata = (upstream_data*) current->user;
printf("ID = %u loc = %u\n",current->ID,i);
for(j=0;j<current->numparents;j++)
{
printf("upstream: %u\n",updata->num_upstream[j]);
for(l=0;l<updata->num_upstream[j];l++)
	printf("%u ",updata->upstream[j][l]);
printf("\n");
}
printf("+++++++\n");
}
*/
		//Clean up
		free(influenced_gauges);
		free(distance);
	}
}


//Deletes any unneeded upstream link information. Use this after system partitioning is determined.
void Clean_Upstream_Links(asynchsolver* asynch)
{
	Link** sys = asynch->sys;
	unsigned int N = asynch->N,i,j;
	int* assignments = asynch->assignments;
	short int* getting = asynch->getting;
	upstream_data* data;

	for(i=0;i<N;i++)
	{
		if(assignments[i] != my_rank && !getting[i])
		{
			data = (upstream_data*) (sys[i]->user);
			if(data)
			{
				for(j=0;j<sys[i]->numparents;j++)
					if(data->upstream[j])	free(data->upstream[j]);
				if(data->fit_states)	free(data->fit_states);
				if(data->fit_to_universal)	free(data->fit_to_universal);
				free(data->upstream);
				free(data->num_upstream);
				free(data);
				sys[i]->user = NULL;
			}
		}
	}
}

void Free_Upstream_Links(asynchsolver* asynch)
{
	upstream_data *data;
	Link** sys = asynch->sys;
	unsigned int N = asynch->N,i,j;

	for(i=0;i<N;i++)
	{
		data = (upstream_data*) (sys[i]->user);
		if(data)
		{
			for(j=0;j<sys[i]->numparents;j++)
				if(data->upstream[j])	free(data->upstream[j]);
			free(data->upstream);
			free(data->num_upstream);
			free(data->fit_states);
			free(data->fit_to_universal);
			free(data);
			sys[i]->user = NULL;
		}
	}
}

//Read into memory the times and discharges stored in a .dat file.
double*** ReadSolution(char filename[],unsigned int** id_to_loc,unsigned int N,unsigned int* numlinks,unsigned int** ids,unsigned int** locs,unsigned int** numsteps)
{
	unsigned int i,j,k,dim;
	double*** data = NULL;

	if(my_rank == 0)
	{
		FILE* file = fopen(filename,"r");
		if(!file)
		{
			printf("Error reading true solution. File %s not found.\n",filename);
			abort();
		}

		//Setup data structures
		fscanf(file,"%u%u",numlinks,&dim);
		*ids = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
		*locs = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
		*numsteps = (unsigned int*) calloc(*numlinks,sizeof(unsigned int));
		data = (double***) malloc(*numlinks*sizeof(double**));

		//Read in the file
		for(i=0;i<*numlinks;i++)
		{
			fscanf( file,"%u %u",&((*ids)[i]),&((*numsteps)[i]) );
			data[i] = (double**) malloc((*numsteps)[i]*sizeof(double*));
			for(j=0;j<(*numsteps)[i];j++)
			{
				data[i][j] = (double*) malloc(2*sizeof(double));
				fscanf(file,"%lf %lf",&(data[i][j][0]),&(data[i][j][1]));	//time and discharge
				for(k=2;k<dim;k++)	fscanf(file,"%*f");
			}
		}

		//Find locations from ids
		for(i=0;i<*numlinks;i++)
			(*locs)[i] = find_link_by_idtoloc((*ids)[i],id_to_loc,N);

		//Cleanup
		fclose(file);
	}

/*
	//Send data to all procs
	MPI_Bcast(numlinks,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	if(my_rank != 0)
	{
		*ids = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
		*locs = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
		*numsteps = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
		data = (double***) malloc(*numlinks*sizeof(double**));
	}
	MPI_Bcast(*ids,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(*locs,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(*numsteps,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	for(i=0;i<*numlinks;i++)
	{
		if(my_rank != 0)
			data[i] = (double**) malloc((*numsteps)[i]*sizeof(double*));
		for(j=0;j<(*numsteps)[i];j++)
		{
			if(my_rank != 0)	data[i][j] = (double*) malloc(2*sizeof(double));
			MPI_Bcast(data[i][j],2,MPI_DOUBLE,0,MPI_COMM_WORLD);
		}
	}
*/

	MPI_Bcast(numlinks,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	if(my_rank != 0)
	{
		*numsteps = NULL;
		*ids = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
		*locs = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
	}
	MPI_Bcast(*ids,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(*locs,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	return data;
}


//Returns the discharges for all links with observations in data at time t.
void FindAllDischarges(double*** data,double t,unsigned int numlinks,unsigned int* numsteps,double* d)
{
	unsigned int i,j,min,max;

	if(my_rank == 0)
	{
		for(j=0;j<numlinks;j++)
		{
			i = numsteps[j]/2;
			max = numsteps[j];
			min = 0;
			if(t > data[j][numsteps[j]-1][0])
				i = max-1;
			else
			{
				while( t*.99999 > data[j][i][0] || t*1.00001 < data[j][i][0] )
				{
					if(data[j][i][0] < t)	min = i+1;
					else			max = i;
					i = (max+min)/2;

					if(min >= max)
					{
						printf("Time %f not found for data %u.\n",t,j);
						break;
					}
				}
			}

			d[j] = data[j][i][1];
		}
	}

	MPI_Bcast(d,numlinks,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

//Gives a list of all locations with a downstream gauge, sorted by location.
//This allocates an appropriate amount of space for above_gauges and bool_above_gauges.
//bool_above_gauges[i] == 1 if sys[i] is above a gauge. 0 if not.
//above_gauges[i] is the list of link locations above a gauge.
//Returns the total number of gauges above a gauge.
unsigned int GaugeDownstream(asynchsolver* asynch,unsigned int** above_gauges,short int** bool_above_gauges,unsigned int* gauges,unsigned int numdata)
{
	if(!above_gauges || !bool_above_gauges)	return 0;
	unsigned int N = asynch->N,num_above,i,j;
	Link *current,*next,**sys = asynch->sys;
	(*bool_above_gauges) = (short int*) calloc(N,sizeof(short int));

	if(my_rank == 0)
	{
		//Set the location with gauges
		for(i=0;i<numdata;i++)
			(*bool_above_gauges)[gauges[i]] = 1;
		num_above = numdata;

		//Trickle down from each external link, until a gauge is reached, marking all links
		for(i=0;i<N;i++)
		{
			if(sys[i]->numparents == 0)
			{
				current = sys[i];
				next = current;

				//See if there is a gauge
				while(next && (*bool_above_gauges)[next->location] == 0)	next = next->c;

				if(next)	//Gauge found
				{
					for(;current != next;current = current->c)
					{
						(*bool_above_gauges)[current->location] = 1;
						num_above++;
					}
				}
			}
		}

		//Setup above_gauges
		*above_gauges = (unsigned int*) malloc(num_above*sizeof(unsigned int));
		j = 0;
		for(i=0;i<N;i++)
		{
			if((*bool_above_gauges)[i])
				(*above_gauges)[j++] = i;
		}
	}

	MPI_Bcast(*bool_above_gauges,N,MPI_SHORT,0,MPI_COMM_WORLD);
	MPI_Bcast(&num_above,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	if(my_rank != 0)	*above_gauges = (unsigned int*) malloc(num_above*sizeof(unsigned int));
	MPI_Bcast(*above_gauges,num_above,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	return num_above;
}

//Scales discharge upstream by area from the gauges
//!!!! Going to need areas when calling this with multiple procs !!!!
int AdjustDischarges_Scale(asynchsolver* asynch,unsigned int* data_locs,double* d,unsigned int numdata,double* x,unsigned int allstates,unsigned int problem_dim)
{
	//Unpack
	UnivVars* GlobalVars = asynch->GlobalVars;
	unsigned int N = asynch->N,**id_to_loc = asynch->id_to_loc;
	int *assignments = asynch->assignments;
	Link** sys = asynch->sys;
	unsigned int area_idx = GlobalVars->area_idx;

	unsigned int i,j,loc,num_upstream,prev_loc,curr_loc;
	unsigned int *upstream = (unsigned int*) malloc(N*sizeof(unsigned int));
	short int *locs_set = (short int*) calloc(N,sizeof(short int));	//1 if the discharge is to be changed, 0 if not
	double *locs_newq = (double*) malloc(N*sizeof(double));	//The new value for the discharge at location i
	unsigned int *counter = (unsigned int*) malloc(N*sizeof(unsigned int));
	double* upstream_areas = (double*) malloc(N*sizeof(double));	//The upstream area of the link location i
	Link* current;
	double ratio;

	//Kinda crappy, but get the upstream area for each link
	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank)	upstream_areas[i] = sys[i]->params.ve[area_idx];
		MPI_Bcast(&(upstream_areas[i]),1,MPI_DOUBLE,assignments[i],MPI_COMM_WORLD);
	}

	//Initialize the new discharges to a negative value
	for(i=0;i<N;i++)	locs_newq[i] = -1.0;

	//Store gauge readings
	for(i=0;i<numdata;i++)
	{
		locs_set[data_locs[i]] = 1;
		locs_newq[data_locs[i]] = d[i];
	}

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
				ratio = locs_newq[current->location] / upstream_areas[current->location];
				//ratio = locs_newq[current->location] / current->params->ve[area_idx];
				for(j=0;j<num_upstream;j++)
				{
					locs_set[upstream[j]] = 1;
					locs_newq[upstream[j]] = upstream_areas[upstream[j]] * ratio;
					//locs_newq[upstream[j]] = sys[upstream[j]]->params->ve[area_idx] * ratio;
				}
			}
		}
	}
/*
	//Set the discharges downstream from each gauge. This is for locations with no downstream gauages.
	//This follows each gauge downstream until a link is found with locs_set or an outlet.
	//counter is set to help track the closest gauge.
	for(i=0;i<N;i++)	counter[i] = N;

	for(i=0;i<numdata;i++)
	{
		loc = data_locs[i];
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
					locs_newq[curr_loc] = upstream_areas[curr_loc] * locs_newq[loc] / upstream_areas[loc];
					//locs_newq[curr_loc] = current->params->ve[area_idx] * locs_newq[loc] / sys[loc]->params->ve[area_idx];
				}
				prev_loc = curr_loc;
				current = current->c;
				if(current)	curr_loc = current->location;
				else		break;
			}
		}
	}

	for(i=0;i<N;i++)
		if(counter[i] < N)	locs_set[i] = 1;
*/


	//Set the determined discharge. If a link's discharge was not determined by the above process, then it lies in a totally ungauged basin.
	for(i=0;i<N;i++)
	{
		if(locs_set[i])
			x[problem_dim*i] = locs_newq[i];
	}

	//Clean up
	free(locs_set);
	free(upstream);
	free(locs_newq);
	free(counter);
	free(upstream_areas);
	return 0;
}

//Reads a .das file
AssimData* Init_AssimData(char* assim_filename,asynchsolver* asynch)
{
	unsigned int **id_to_loc = asynch->id_to_loc,N = asynch->N,string_size = asynch->GlobalVars->string_size,i,id,n,dropped;
	int errorcode,valsread;
	AssimData* Assim;
	FILE* inputfile;
	char end_char;
	unsigned int buff_size = string_size + 20;
	char* linebuffer = (char*) malloc(buff_size*sizeof(char));
	MPI_Barrier(MPI_COMM_WORLD);

	if(!id_to_loc)
	{
		if(my_rank == 0)	printf("Error: Topology data must be processed before starting data assimilation routines.\n");
		return NULL;
	}

	if(my_rank == 0)
	{
		//Open file
		inputfile = fopen(assim_filename,"r");
		errorcode = 0;
		if(!inputfile)
		{
			printf("[%i]: Error opening assimilation file %s.\n",my_rank,assim_filename);
			errorcode = 1;
		}
	}

	//Check if the assimilation file was openned
	MPI_Bcast(&errorcode,1,MPI_INT,0,MPI_COMM_WORLD);
	if(errorcode)	return NULL;

	//Create assim
	Assim = (AssimData*) malloc(sizeof(AssimData));
	Assim->conninfo = NULL;
	Assim->data_locs = NULL;
	Assim->id_to_assim = NULL;
	Assim->db_filename = (char*) malloc(string_size*sizeof(char));

	//Read the .dbc file
	ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%s",Assim->db_filename);
	if(ReadLineError(valsread,1,"assimilation dbc filename"))	return NULL;
	Assim->conninfo = ReadDBC(Assim->db_filename,string_size);
	if(!Assim->conninfo)	return NULL;

	//Read in rest of data
	ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%u",&(Assim->steps_to_use));
	if(ReadLineError(valsread,1,"number of observations to use"))	return NULL;

	ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%lf",&(Assim->inc));
	if(ReadLineError(valsread,1,"time resolution of observations"))	return NULL;

	//ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
	//valsread = sscanf(linebuffer,"%lf",&(Assim->forecast_window));
	//if(ReadLineError(valsread,1,"forecast window"))	return NULL;

	ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%u",&(Assim->least_squares_iters));
	if(ReadLineError(valsread,1,"least squares iterations"))	return NULL;

	//Read ending mark
	ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%c",&end_char);
	if(ReadLineError(valsread,1,"ending mark"))	return NULL;

	//Clean up
	if(my_rank == 0)
		fclose(inputfile);
	free(linebuffer);
	MPI_Barrier(MPI_COMM_WORLD);
	if(end_char != '#')
	{
		if(my_rank == 0)	printf("Error: Ending mark not seen in %s.\n",assim_filename);
		return NULL;
	}

	return Assim;
}

//Trashes an AssimData object
void Free_AssimData(AssimData** assim)
{
	unsigned int i;
	if((*assim)->id_to_assim)
	{
		for(i=0;i<(*assim)->numdata;i++)	free((*assim)->id_to_assim[i]);
		free((*assim)->id_to_assim);
	}
	if((*assim)->conninfo)	ConnData_Free((*assim)->conninfo);
	if((*assim)->data_locs)	free((*assim)->data_locs);
	if((*assim)->db_filename)	free((*assim)->db_filename);
	free(*assim);
	*assim = NULL;
}

int Download_Gauge_IDs(asynchsolver* asynch,AssimData* Assim)
{
	int errorcode = 0;
	unsigned int i,n,dropped,*gauged_ids = NULL,**id_to_loc = asynch->id_to_loc,N = asynch->N;
	char query[ASYNCH_MAX_QUERY_LENGTH];
	PGresult *res;

	//Get link ids with gauges
	if(my_rank == 0)
	{
		ConnectPGDB(Assim->conninfo);
		sprintf(query,Assim->conninfo->queries[0]);
		res = PQexec(Assim->conninfo->conn,query);
		if(CheckResError(res,"getting list of gauge link ids"))
			errorcode = 1;
		else
		{
			dropped = 0;
			n = PQntuples(res);
			Assim->numdata = n;
			Assim->data_locs = (unsigned int*) malloc(n*sizeof(unsigned int));
			gauged_ids = (unsigned int*) malloc(n*sizeof(unsigned int));
			for(i=0;i<n;i++)
			{
				gauged_ids[i-dropped] = atoi(PQgetvalue(res,i,0));
				Assim->data_locs[i-dropped] = find_link_by_idtoloc(gauged_ids[i],id_to_loc,N);
				if(Assim->data_locs[i-dropped] > N)
				{
					printf("Warning: Ignoring gauge at link id %u. No link with this id is in the network.\n",gauged_ids[i-dropped]);
					dropped++;
				}
			}

			//Resize the data_locs array
			Assim->numdata -= dropped;
			Assim->data_locs = (unsigned int*) realloc(Assim->data_locs,Assim->numdata * sizeof(unsigned int));
		}

		PQclear(res);
		DisconnectPGDB(Assim->conninfo);
	}

	//Give the locations to the other procs
	MPI_Bcast(&errorcode,1,MPI_INT,0,MPI_COMM_WORLD);
	if(!errorcode)
	{
		MPI_Bcast(&(Assim->numdata),1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		n = Assim->numdata;
		if(my_rank != 0)
		{
			Assim->data_locs = (unsigned int*) malloc(n*sizeof(unsigned int));
			gauged_ids = (unsigned int*) malloc(n*sizeof(unsigned int));
		}
		MPI_Bcast(Assim->data_locs,n,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		MPI_Bcast(gauged_ids,n,MPI_UNSIGNED,0,MPI_COMM_WORLD);

		//Build id_to_assim
		Assim->id_to_assim = (unsigned int**) malloc(n*sizeof(unsigned int*));	//!!!! Should really have an array (should REALLY be a hash table) !!!!
		for(i=0;i<n;i++)	Assim->id_to_assim[i] = (unsigned int*) malloc(2*sizeof(unsigned int*));
		for(i=0;i<n;i++)
		{
			Assim->id_to_assim[i][0] = gauged_ids[i];
			Assim->id_to_assim[i][1] = i;
		}
		merge_sort_ids(Assim->id_to_assim,n);
	}

	if(gauged_ids)	free(gauged_ids);
	return errorcode;
}

//Download gauge data and store them in d. d orders the data by time, then location.
//Data are obtained from times starting at background_time_unix and ending at steps_to_use*inc+background_time_unix.
//Returns 0 if everything went as planned. Returns 1 if some data was not available. Returns -1 if an error occurred.
int GetObservationsDB(AssimData* Assim,unsigned int **id_loc_loc,unsigned int N,unsigned int background_time_unix,double* d)
{
	unsigned int i,j,n,end_time_unix,*data_locs = Assim->data_locs,idx,inc_secs = (unsigned int)(Assim->inc*60.0 + 1e-3),unix_t,id,numdata = Assim->numdata;
	unsigned int steps_to_use = Assim->steps_to_use,current_time;
	int errorcode = 0;
	PGresult *res;
	char *query = Assim->conninfo->query;

	//Reset d
	for(i=0;i<steps_to_use*numdata;i++)	d[i] = -1.0;

	//Download the data
	if(my_rank == 0)
	{
		if(ConnectPGDB(Assim->conninfo))	errorcode = -1;
		else
		{
			end_time_unix = background_time_unix + steps_to_use * 60 * (int) (Assim->inc + 1e-3);
			sprintf(query,Assim->conninfo->queries[1],background_time_unix,end_time_unix);	//Assumes data sorted by time, then link id (!!!! No, is it the reverse? !!!!)
			res = PQexec(Assim->conninfo->conn,query);
			if(CheckResError(res,"downloading gauge data") == 0)
			{
				//Unpack the data
				n = PQntuples(res);
				for(i=0;i<n;)
				{
					id = atoi(PQgetvalue(res,i,0));
/*
printf("Unpacking data for id %u...\n",id);
int stopper = 1;
while(id == 204046 && stopper == 1)
{
sleep(5);
}
*/
					idx = find_link_by_idtoloc(id,Assim->id_to_assim,numdata);	//!!!! What if no data for a link is available? Doing nothing keeps the values the same. Or zero... !!!!
					if(idx < numdata)
					{
						while(i<n && atoi(PQgetvalue(res,i,0)) == id)
						{
							//Check if the first time needed is available
							unix_t = atoi(PQgetvalue(res,i,1));
							if(background_time_unix == unix_t)
							{
								//d[idx*steps_to_use] = atof(PQgetvalue(res,i,2));
								d[idx] = atof(PQgetvalue(res,i,2));
								i++;
							}
							else	//For now, just take from the next available reading
							{
								d[idx] = atof(PQgetvalue(res,i,2));
								//printf("!!!! Error: didn't get a gauge value for id = %u at time %u. !!!!\n",id,background_time_unix);
								//MPI_Abort(MPI_COMM_WORLD,1);
							}

							//i++;
							current_time = 1;
							unix_t = background_time_unix + inc_secs;
							//unix_t += inc_secs;
							//for(;i<n && atoi(PQgetvalue(res,i,0)) == id;i++)
							while(i<n && atoi(PQgetvalue(res,i,0)) == id)
							{
								if(atoi(PQgetvalue(res,i,1)) == unix_t)
								{
									//d[idx*steps_to_use+current_time] = atof(PQgetvalue(res,i,2));
									d[idx+current_time*numdata] = atof(PQgetvalue(res,i,2));
									current_time++;
									unix_t += inc_secs;
									i++;
								}
								else if(atoi(PQgetvalue(res,i,1)) > unix_t)	//Not available!
								{
									//Just use the last observation
									//d[idx*steps_to_use+current_time] = atof(PQgetvalue(res,i-1,2));
									//d[idx+current_time*numdata] = atof(PQgetvalue(res,i-1,2));
									d[idx+current_time*numdata] = d[idx+(current_time-1)*numdata];
									current_time++;
									unix_t += inc_secs;
									errorcode = 1;
								}
								else	//The time step for this location does not match time_inc
								{
									i++;	//This does NOT use the closest gauge reading. It basically ignores all off-timestep observations
								}
							}
						}
					}
					else	//ID not a gauged location. So skip it.
					{
						i++;
						while(i<n && atoi(PQgetvalue(res,i,0)) == id)	i++;
					}
				}
			}
			else
				errorcode = -1;

			//Clean up
			PQclear(res);
			DisconnectPGDB(Assim->conninfo);
		}
	}

	//Send the data to everyone
	MPI_Bcast(&errorcode,1,MPI_INT,0,MPI_COMM_WORLD);
	if(errorcode > -1)
		MPI_Bcast(d,steps_to_use*numdata,MPI_DOUBLE,0,MPI_COMM_WORLD);
if(my_rank == 0)
{
printf("At time %u, got:\n",background_time_unix);
for(i=0;i<steps_to_use*numdata;i++)
	printf("%f ",d[i]);
printf("\n+++++++++++++++++++++++\n");
}
	return errorcode;
}


//Checks which dicharge values are having problems converging. Everything upstream is cut in half.
//!!!! This is a little sloppy with overlapping domains of influence. Could reduce by more than 1/2. Also broadcasting static information. !!!!
int ReduceBadDischargeValues(Link** sys,int* assignments,unsigned int N,double* d_full,double* q,unsigned int steps_to_use,unsigned int* data_locs,unsigned int numdata,double* x_start,unsigned int assim_dim,double limit)
{
	unsigned int i,j,k,loc,num_links;
	int value = 0;
	Link* current;
	double new_diff,old_diff,factor = 0.8;
	upstream_data* updata;

	for(i=0;i<numdata;i++)
	{
		if(d_full[i] < 0.0)
		{
			if(my_rank == 0)
				printf("!!!! Skipping scaling of link %u. Observation is %f. !!!!\n",data_locs[i],d_full[i]);
			continue;
		}

		current = sys[data_locs[i]];

		//Check how far off each gauge is
		//old_diff = d_full[i] - q[i];
		old_diff = q[i] - d_full[i];
		for(j=1;j<steps_to_use;j++)
		{
			//new_diff = d_full[i+j*numdata] - q[i+j*numdata];
			new_diff = q[i+j*numdata] - d_full[i+j*numdata];
			if(old_diff < new_diff)	//Getting worse...
				old_diff = new_diff;
			else
				break;
		}

		//if(j == steps_to_use && old_diff > limit)	//Got progressively worse, and ended badly. Make a change!
		if(j == steps_to_use)
		{
			value = 1;

			//Calculate the scaling factor
			factor = d_full[i+(steps_to_use-1)*numdata] / q[i+(steps_to_use-1)*numdata];

			//Reduce q at the gauge
			x_start[data_locs[i] * assim_dim] *= factor;

			if(assignments[data_locs[i]] == my_rank)
			{
				printf("[%i] !!!! Reducing discharges above link %u. Factor is %f. !!!!\n",my_rank,current->ID,factor);
				updata = (upstream_data*) current->user;

				//Sum up the total number of links upstream from current
				num_links = 0;
				for(j=0;j<current->numparents;j++)	num_links += updata->num_upstream[j];
				MPI_Bcast(&num_links,1,MPI_UNSIGNED,my_rank,MPI_COMM_WORLD);

				for(j=0;j<current->numparents;j++)
				{
					for(k=0;k<updata->num_upstream[j];k++)
					{
						MPI_Bcast(&(updata->upstream[j][k]),1,MPI_UNSIGNED,my_rank,MPI_COMM_WORLD);
						x_start[updata->upstream[j][k] * assim_dim] *= factor;
					}
				}
			}
			else
			{
				MPI_Bcast(&num_links,1,MPI_UNSIGNED,assignments[data_locs[i]],MPI_COMM_WORLD);

				for(j=0;j<num_links;j++)
				{
					MPI_Bcast(&loc,1,MPI_UNSIGNED,assignments[data_locs[i]],MPI_COMM_WORLD);
					x_start[loc * assim_dim] *= factor;
				}
			}
		}
	}

	return value;
}


//This creates a snapshot of the model states only (i.e. no variational equation states are saved)
int SnapShot_ModelStates(asynchsolver* asynch,unsigned int problem_dim)
{
	Link** sys = asynch->sys;
	UnivVars* GlobalVars = asynch->GlobalVars;
	int *assignments = asynch->assignments;
	unsigned int i,j,N = asynch->N;
	FILE* output;
	double buffer[ASYNCH_MAX_DIM];

	if(my_rank == 0)	//Creating the file
	{
		output = fopen(GlobalVars->dump_loc_filename,"w");
		if(output == NULL)
		{
			printf("[%i]: Error opening file %s.\n",my_rank,GlobalVars->dump_loc_filename);
			i = 1;
		}
		else	i = 0;
		MPI_Bcast(&i,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		if(i)	return 1;

		fprintf(output,"%hu\n%u\n0.0\n\n",GlobalVars->type,N);

		for(i=0;i<N;i++)
		{
			if(assignments[i] != 0)
				MPI_Recv(buffer,problem_dim,MPI_DOUBLE,assignments[i],i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			else
				for(j=0;j<problem_dim;j++)	buffer[j] = sys[i]->list->tail->y_approx.ve[j];

			fprintf(output,"%u\n",sys[i]->ID);
			for(j=0;j<problem_dim;j++)	fprintf(output,"%.6e ",buffer[j]);
			fprintf(output,"\n");
		}

		fclose(output);
	}
	else			//Sending data to proc 0
	{
		//Check for error
		MPI_Bcast(&i,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		if(i)	return 1;

		for(i=0;i<N;i++)
		{
			if(assignments[i] == my_rank)
			{
				for(j=0;j<problem_dim;j++)	buffer[j] = sys[i]->list->tail->y_approx.ve[j];
				MPI_Send(buffer,problem_dim,MPI_DOUBLE,0,i,MPI_COMM_WORLD);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return 0;
}

int GaugeDataAvailable(AssimData* Assim,unsigned int start_time,unsigned int end_time)
{
	int isnull = 1,numdata,numvalues,expect;
	PGresult* res;
	char* query = NULL;

	if(my_rank == 0)
	{
		query = (char*) malloc(1024*sizeof(char));
		ConnectPGDB(Assim->conninfo);
		numdata = Assim->numdata;

		//Get number of observations available
		sprintf(query,Assim->conninfo->queries[2],start_time,end_time);
		res = PQexec(Assim->conninfo->conn,query);
		if(CheckResError(res,"checking for new observation data"))	goto error;
		numvalues = atoi(PQgetvalue(res,0,0));
		PQclear(res);

		//Is a good number of observations available? Waiting for 80%.
		expect = (Assim->steps_to_use * numdata * 8)/10;
		if(numvalues >= expect)	isnull = 0;

		printf("%i/%i new observations available. Expecting at least %i.\n",numvalues,Assim->steps_to_use * numdata,expect);

		error:
		DisconnectPGDB(Assim->conninfo);
	}

	MPI_Bcast(&isnull,1,MPI_INT,0,MPI_COMM_WORLD);
	if(query)	free(query);
	return isnull;
}

