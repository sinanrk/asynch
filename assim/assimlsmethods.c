#include "assimlsmethods.h"

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
			for(j=0;j<problem_dim;j++)	current->list->head->y_approx->ve[j] = backup[i*problem_dim+j];
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



//Read into memory the times and discharges stored in a .dat file.
double*** ReadSolution(char filename[],unsigned int** id_to_loc,unsigned int N,unsigned int* numlinks,unsigned int** ids,unsigned int** locs,unsigned int** numsteps)
{
	unsigned int i,j,k,dim;
	double*** data;

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

/*
//!!!! Add noise !!!!
for(i=0;i<*numlinks;i++)
{
	for(j=0;j<(*numsteps)[i];j++)
	{
		data[i][j][1] += ((rand() % 2000) - 1000)/1000.0;
		if(data[i][j][1] < 1e-6)	data[i][j][1] = 1e-6;
	}
}
*/

		//Cleanup
		fclose(file);
	}

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

	return data;
}


//Returns the discharges for all links with observations in data at time t.
void FindAllDischarges(double*** data,double t,unsigned int numlinks,unsigned int* numsteps,double* d)
{
	unsigned int i,j,min,max;

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
				else				max = i;
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

/*
//Returns the discharge for link id at time t.
//!!!! Could this be given the location instead of id? !!!!
double FindDischarge(double*** data,unsigned int id,double t,unsigned int numlinks,unsigned int* ids,unsigned int* numsteps)
{
	unsigned int i,min,max,found_loc;

	//Find the location of the id
	for(i=0;i<numlinks;i++)
	{
		if(ids[i] == id)
		{
			found_loc = i;
			break;
		}
	}
	if(i == numlinks)
	{
		printf("Link id %u not found.\n",id);
		return -1.0;
	}

	//Find the time
	i = numsteps[found_loc]/2;
	max = numsteps[found_loc];
	min = 0;
	if(t > data[found_loc][numsteps[found_loc]-1][0])
		i = max-1;
	else
	{
		while( t*.99999 > data[found_loc][i][0] || t*1.00001 < data[found_loc][i][0] )
		{
			if(data[found_loc][i][0] < t)	min = i;
			else				max = i;
			i = (max+min)/2;

			if(min == max)
			{
				printf("Time %f not found.\n",t);
				return -1.0;
			}
		}
	}

	return data[found_loc][i][1];
}
*/


