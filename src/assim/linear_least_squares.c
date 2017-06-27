#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <mpi.h>

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#if defined(HAVE_PETSC)
#include <petsc.h>
#endif

#if !defined(_MSC_VER)
#define ASYNCH_SLEEP(seconds) sleep(seconds)
#else
#include <windows.h>
#define ASYNCH_SLEEP(seconds) Sleep((seconds) * 1000)
#endif

#include <asynch_interface.h>

//Internal stuffs
#include <blas.h>
#include <comm.h>
#include <system.h>

#include <assim/models.h>
#include <assim/structs.h>
#include <assim/linear_least_squares.h>


/// This computes the least squares fit assuming the background and analysis difference is linear in the innovations. 
/// HM is (num_obs*steps_to_use) X allstates_needed
/// HM_els is 1 X allstates (i.e. a 1D array)
/// HM_buffer is 1 X allstates_needed
/// \param asynch The asynch solver instance
/// \param asynch The assimilation workspace
int SolveSysLS(AsynchSolver* asynch, AssimWorkspace* ws, double* q)
{
    unsigned int i, j,/*k,m,*/n/*,l,counter*/;

    //Unpack ptr
    unsigned int N = asynch->N;
    Link* sys = asynch->sys;
    GlobalVars* globals = asynch->globals;
    int* assignments = asynch->assignments;
    short int* getting = asynch->getting;
    //Mat *HM = ws->HM, *HTH = ws->HTH, *HMTR = ws->HMTR;
    //Vec *RHS = ws->RHS, d, *x = ws->x, *B = ws->B, *R = ws->R;
    //KSP *ksp = ws->ksp;
    Vec d;
    unsigned int *obs_locs = ws->obs_locs, assim_dim = ws->assim_dim;
    unsigned int problem_dim = ws->problem_dim, allstates = ws->allstates;
    //unsigned int steps_to_use = ws->steps_to_use, num_obs = ws->num_obs;
    int *HM_col_indices = ws->HM_col_indices, *d_indices = ws->d_indices;
    double t_b = ws->t_b;
    unsigned int allstates_needed = ws->allstates_needed;
    double /**RHS_els,*/*x_start = ws->x_start, *HM_buffer = ws->HM_buffer;
    //unsigned int max_or_steps = steps_to_use;
    AsynchModel* custom_model = asynch->model;
    unsigned int *vareq_shift = ws->vareq_shift, *inv_vareq_shift = ws->inv_vareq_shift;

    double start = MPI_Wtime();

    unsigned int num_total_obs = ws->num_steps * ws->num_obs;

    //Build a vector structure for d !!!! Obviously, this needs to not happen... !!!!
    //if(max_or_steps > allstates_needed)	printf("[%i]: Error: max_or_steps > allstates_needed (%u > %u)\n",my_rank,max_or_steps,allstates_needed);
    VecCreateSeq(MPI_COMM_SELF, num_total_obs, &d);	//!!!! This needs to be fixed !!!!
    VecSet(d, 0.0);
    VecSetValues(d, num_total_obs, d_indices, ws->d_full, INSERT_VALUES);
    VecAssemblyBegin(d);
    VecAssemblyEnd(d);

    /*
    //!!!! Try to scale the init conditions !!!!
    AdjustDischarges(asynch,obs_locs,d_els,num_obs,x_start,allstates);
    printf("Adjusted x_start to\n");
    for(i=0;i<allstates;i++)	printf("%.15e ",x_start[i]);
    printf("\n");
    getchar();
    */
    //Initialize the system
    ResetSysLS(sys, N, globals, t_b, x_start, assim_dim, globals->num_forcings, asynch->my_data);

    for (i = 0; i < N; i++)				//!!!! Put this into ResetSysLS? !!!!
        if (assignments[i] == asynch->my_rank || getting[i])
            custom_model->initialize_eqs(
                globals->global_params, globals->num_global_params,
                sys[i].params, globals->num_params,
                sys[i].my->list.head->y_approx, sys[i].dim,
                sys[i].user); //!!!! Should all states be reset? !!!!
            //ReadInitData(globals->global_params,sys[i].params,NULL,0,sys[i].list->head->y_approx,globals->type,sys[i].diff_start,sys[i].no_ini_start,sys[i].user,NULL);	//!!!! Very inefficient. Too many checks. !!!!

    ////Initialize the variational equations
    //unsigned int stack_size = 0;
    //Link **stack = (Link**)calloc(N, sizeof(Link*));
    //for (i = 0; i < N; i++)
    //    if (sys[i].num_parents == 0)
    //        stack[stack_size++] = &sys[i];

    //// Visit from source to outlet
    //unsigned int *visits = (unsigned int *)calloc(N, sizeof(unsigned int));
    //while (stack_size > 0)
    //{
    //    Link *current = stack[stack_size - 1];
    //    UpstreamData* updata = (UpstreamData*)current->user;

    //    // Pop from the stack
    //    stack_size--;

    //    // Increment visit counter of child
    //    if (current->child)
    //        visits[current->child->location]++;

    //    double lambda_1 = globals->global_params.ve[1];
    //    double k_3 = globals->global_params.ve[4];	//[1/min]
    //    double h_b = globals->global_params.ve[6];	//[m]
    //    double S_L = globals->global_params.ve[7];	//[m]
    //    double A = globals->global_params.ve[8];
    //    double B = globals->global_params.ve[9];
    //    double exponent = globals->global_params.ve[10];

    //    double L = current->params.ve[1];	//[m]
    //    double A_h = current->params.ve[2];	//[m^2]
    //    double invtau = current->params.ve[3];	//[1/min]
    //    double k_2 = current->params.ve[4];	//[1/min]
    //    double k_i = current->params.ve[5];	//[1/min]
    //    double c_1 = current->params.ve[6];
    //    double c_2 = current->params.ve[7];

    //    VEC y_0 = current->list->head->y_approx;
    //    double q = y_0.ve[0];		//[m^3/s]
    //    double s_p = y_0.ve[1];	//[m]
    //    double s_t = y_0.ve[2];	//[m]
    //    double s_s = y_0.ve[3];	//[m]

    //    //unsigned int i;
    //    unsigned int offset = 4;

    //    y_0.ve[offset] = 1.;  //dq/dq_0

    //    //A few calculations...
    //    double q_to_lambda_1 = pow(q, lambda_1);
    //    double q_to_lambda_1_m1 = pow(q, lambda_1 - 1.0);// (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);

    //    //Discharge
    //    double inflow = 0.0;
    //    for (unsigned int i = 0; i < updata->num_parents; i++)
    //        inflow += updata->parents[i]->list->head->y_approx.ve[0];

    //    //Compute partial derivatives (local variables)
    //    double dfq_dq = lambda_1 * invtau * q_to_lambda_1_m1 * (-q + c_2*(k_2*s_p + k_3*s_s) + inflow) - invtau * q_to_lambda_1;

    //    //Compute partial derivatives (upstreams variables)
    //    double dfq_dupq = invtau*q_to_lambda_1;

    //    //y_0.ve[offset] = dfq_dq * y_0.ve[offset]; //dq/dq_0

    //    unsigned int j = 0, p = 0;
    //    // For every upstream links
    //    for (unsigned int i = 0, j = 0; i < updata->num_upstreams; i++, j++)
    //    {
    //        unsigned int np = p + 1;

    //        // If switch to next parent
    //        if (np < updata->num_parents && updata->upstreams[i] == updata->parents[np])
    //        {
    //            p++;
    //            j = 0;
    //        }

    //        unsigned int current_idx = offset + i + 1;
    //        unsigned int parent_idx = offset + j;
    //        VEC y_p = updata->parents[p]->list->head->y_approx;

    //        assert(current_idx < y_0.dim);
    //        assert(parent_idx < y_p.dim);
    //        y_0.ve[current_idx] = dfq_dupq * y_p.ve[parent_idx] + dfq_dq * y_0.ve[current_idx]; //q, upq
    //    }

    //    if (current->child && visits[current->child->location] == current->child->num_parents)
    //        stack[stack_size++] = current->child;
    //}

    for (i = 0; i < asynch->globals->num_forcings; i++)
    {
        if (asynch->forcings[i].flag == 3)	//!!!! Recurring and binary files need this too !!!!
        {
            //printf("Setting to %u and %u\n",asynch->forcings[i]->first_file,asynch->forcings[i]->last_file);
            Asynch_Set_Forcing_State(asynch, i, t_b, asynch->forcings[i].first_file, asynch->forcings[i].last_file);
        }
    }

    //Advance the system and extract the HM matrix
    for (i = 0; i < ws->num_steps; i++)	//!!!! Start at i=1? For i = 0, I don't think we need to set anything... !!!!		//HM here holds the values of M that are needed
    {
        globals->t = 0.0;

        if (i > 0)
        {
            // Adjust the end of the simulation
            globals->maxtime = t_b + i * ws->obs_time_step;

            MPI_Barrier(MPI_COMM_WORLD);
            double start = MPI_Wtime();

            Asynch_Advance(asynch, 0);

            MPI_Barrier(MPI_COMM_WORLD);
            double stop = MPI_Wtime();

            if (asynch->my_rank == 0)
                printf("Time for advance to time %f: %.0f\n", globals->maxtime, stop - start);
        }


        //printf("ID = %u, t = %e\n",sys[obs_locs[0]]->ID,sys[obs_locs[0]]->last_t);
        //Print_Vector(sys[obs_locs[0]]->list->tail->y_approx);
        //printf("**********\n");

        //Build HM
        for (j = 0; j < ws->num_obs; j++)
        {
            Link *current = &sys[obs_locs[j]];	//!!!! Assumes only discharges !!!!
            int owner = assignments[obs_locs[j]];
            bool is_my_link = (owner == asynch->my_rank);

            UpstreamData *updata = (UpstreamData*)(current->user);
            memset(HM_buffer, 0, allstates_needed * sizeof(double));

            //From my link
            if (is_my_link)
            {
                //Pull out needed data
                for (n = 0; n < updata->num_fit_states; n++)
                {
                    if (asynch->verbose)
                        printf("ID = %u | Loading %e (from %u) into spot %u\n",
                            current->ID,
                            current->my->list.tail->y_approx[updata->fit_states[n]],
                            updata->fit_states[n],
                            vareq_shift[updata->fit_to_universal[n]]);

                    assert(updata->fit_states[n] < current->dim);
                    HM_buffer[vareq_shift[updata->fit_to_universal[n]]] = current->my->list.tail->y_approx[updata->fit_states[n]];
                }

                //Extract calculationed q's (Just needed for testing. Maybe...)
                q[i * ws->num_obs + j] = current->my->list.tail->y_approx[0];
            }

            //MPI_Bcast(HM_buffer, allstates_needed, MPI_DOUBLE, owner, MPI_COMM_WORLD);	//!!!! Only proc 0 needs this !!!!
            if (asynch->my_rank == 0)
            {
                MPI_Reduce(MPI_IN_PLACE, HM_buffer, allstates_needed, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

#if !defined(NDEBUG)
                unsigned int k;
                for (k = 0; k < allstates_needed; k++)
                    if (HM_buffer[k] != 0.)
                        break;

                assert(k < allstates_needed);
#endif

                unsigned int row_idx = i * ws->num_obs + j;
                MatSetValues(ws->HM, 1, &row_idx, allstates_needed, HM_col_indices, HM_buffer, INSERT_VALUES);
            }
            else
                MPI_Reduce(HM_buffer, NULL, allstates_needed, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



            MPI_Bcast(&(q[i * ws->num_obs + j]), 1, MPI_DOUBLE, owner, MPI_COMM_WORLD);

            //printf("Got %u, %u\n",allstates_needed,updata->num_fit_states);
            //for(n=0;n<allstates_needed;n++)
            //	printf("%e ",HM_buffer[n]);
            //printf("\n");
            //char r = getchar();                
        }
    }

    double stop = MPI_Wtime();

    if (asynch->my_rank == 0)
        printf("Time for advance to time %f: %.0f\n", globals->maxtime, stop - start);


    /*
        //Zero out any unused rows
        if(i < steps_to_use)
        {
    printf("Starting at %u, going to %u\n",i*num_obs,steps_to_use*num_obs);
            for(j=0;j<allstates;j++)
                HM_els[j] = 0.0;
            for(i=i*num_obs;i<steps_to_use*num_obs;i++)
                MatSetValues(ws->HM,1,&i,allstates,cols_allstates_needed,HM_els,INSERT_VALUES);
        }
    */

    if (asynch->my_rank == 0)
    {
        //Assemble the HM matrix
        MatAssemblyBegin(ws->HM, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(ws->HM, MAT_FINAL_ASSEMBLY);

        if (asynch->verbose)
        {
            printf("Matrix HM\n");
            MatView(ws->HM, PETSC_VIEWER_STDOUT_SELF);
        }

        start = MPI_Wtime();

        //Calculate innovations
        double *buffer = NULL;
        VecGetArray(d, &buffer);
        for (i = 0; i < num_total_obs; i++)
            buffer[i] = buffer[i] - q[i];
        VecRestoreArray(d, &buffer);

        //Build the linear system \f$ A x = rhs \f$
        //HMTR is allstates_needed x (num_obs*max_or_steps)
        //HM is (num_obs*max_or_steps) x allstates_needed

        /// \f$ HMTR = H(y_0)^T R \f$
        /// HMTR is a temporary variable used for rhs and A computation
        MatTranspose(ws->HM, MAT_REUSE_MATRIX, &ws->HMTR);
        MatDiagonalScale(ws->HMTR, NULL, ws->R);

        /// \f$ A = B + H(y_0)^T R H(y_0) \f$
        MatMatMult(ws->HMTR, ws->HM, MAT_REUSE_MATRIX, PETSC_DEFAULT, &ws->HTH);
        MatDiagonalSet(ws->HTH, ws->B, ADD_VALUES);
        /*
                MatTransposeMatMult(*HM,*HM,MAT_REUSE_MATRIX,PETSC_DEFAULT,HTH);
                for(i=0;i<allstates_needed;i++)
                {
                    //MatSetValue(*HTH,i,i,1.0,ADD_VALUES);  //!!!! To skip hillslope !!!!
                    //if(i%2)	MatSetValue(*HTH,i,i,1e3,ADD_VALUES);	//Used for s_p, I think...
                    //if(i%2)	MatSetValue(*HTH,i,i,1e2,ADD_VALUES);	//Used for s_t
                    //else	MatSetValue(*HTH,i,i,1.0,ADD_VALUES);
                }
        */

        /// \f$ rhs = H(y_0)^T R \alpha(y_0^b) \f$
        MatMult(ws->HMTR, d, ws->rhs);
        MatAssemblyBegin(ws->HTH, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(ws->HTH, MAT_FINAL_ASSEMBLY);

        if (asynch->verbose)
        {
            printf("Matrix HTH\n");
            MatView(ws->HTH, PETSC_VIEWER_STDOUT_SELF);
        }

        //MPI_Barrier(MPI_COMM_WORLD);
        double stop = MPI_Wtime();
        if (asynch->my_rank == 0)
            printf("Time for matrix computations: %.0f\n", stop - start);

        //Compute analysis
    //MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        /// \f$ x = y_0 - y_0^b \f$
        KSPSetOperators(ws->ksp, ws->HTH, ws->HTH);     //Maybe not actually necessary
        KSPSolve(ws->ksp, ws->rhs, ws->x);
        KSPConvergedReason reason;
        KSPGetConvergedReason(ws->ksp, &reason);

        if (asynch->my_rank == 0)
            printf("Converged reason: %s\n", KSPConvergedReasons[reason]);

        //MPI_Barrier(MPI_COMM_WORLD);
        stop = MPI_Wtime();
        if (asynch->my_rank == 0)
            printf("Time for inversion: %.0f\n", stop - start);

        if (asynch->verbose)
        {
            printf("Solution x\n");
            VecView(ws->x, PETSC_VIEWER_STDOUT_SELF);
        }

        //Copy new solution to x_start
        VecGetArray(ws->x, &buffer);
        //for(i=0;i<num_above;i++)
        for (i = 0; i < allstates_needed; i++)	//!!!! I think this is right... !!!!
        {
            //printf("i = %u %u\n",i,inv_vareq_shift[i]);
            //ASYNCH_SLEEP(1);
            x_start[inv_vareq_shift[i]] += buffer[i];
            //x_start[above_gauges[i]*assim_dim] += x_els[i];	//!!!! To skip hillslope !!!!
            //for(j=0;j<assim_dim;j++)
            //	x_start[above_gauges[i]*assim_dim+j] += x_els[i*assim_dim+j];
        }
        VecRestoreArray(ws->x, &buffer);
    }

    //Send solution to everyone
    MPI_Bcast(x_start, allstates, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (asynch->verbose && asynch->my_rank == 0)
    {
        //unsigned int idxm[num_obs*steps_to_use];
        //double temp_matptr[(num_obs*steps_to_use*allstates_needed > allstates_needed*allstates_needed) ? num_obs*steps_to_use*allstates_needed : allstates_needed*allstates_needed];
        //for(i=0;i<num_obs*steps_to_use;i++)
        //    idxm[i] = i;

        //printf("x_start\n");
        //for(i=0;i<allstates;i++)
        //    printf("%.15e ",x_start[i]);
        //printf("\n");

        double* buffer;
        printf("difference (x)\n");
        VecGetArray(ws->x, &buffer);
        for (i = 0; i < allstates_needed; i++)
            printf("%.2e ", buffer[i]);
        printf("\n");
        VecRestoreArray(ws->x, &buffer);

        printf("d\n");
        VecGetArray(d, &buffer);
        for (i = 0; i < num_total_obs; i++)
            printf("%.2e ", buffer[i]);
        printf("\n");
        VecRestoreArray(d, &buffer);

        //printf("HM\n");
        //MatGetValues(*HM,num_obs*max_or_steps,idxm,allstates_needed,cols_allstates_needed,temp_matptr);
        //for(i=0;i<num_obs*max_or_steps;i++)
        //{
           // for(j=0;j<allstates_needed;j++)
              //  printf("%.15e ",temp_matptr[i*allstates_needed + j]);
           // printf(";\n");
        //}

        //printf("HTH\n");
        //MatGetValues(*HTH,allstates_needed,cols_allstates_needed,allstates_needed,cols_allstates_needed,temp_matptr);
        //for(i=0;i<allstates_needed;i++)
        //{
           // for(j=0;j<allstates_needed;j++)
              //  printf("%.15e ",temp_matptr[i*allstates_needed + j]);
           // printf(";\n");
        //}
    }

    //if (verbose)
    //{
    //    //Get q's produced from analysis (for testing)
    //    if (my_rank == 0)
    //    {
    //        //printf("q before\n");
    //        //Print_VECTOR(q,num_total_obs);
    //    }
    //    double first_diff = ComputeDiff(d_els, q, num_total_obs);

    //    ResetSysLS(sys, N, globals, t_b, x_start, problem_dim, globals->num_forcings, asynch->my_data);
    //    for (i = 0; i < N; i++)
    //        if (assignments[i] == my_rank || getting[i])
    //            //ReadInitData(globals->global_params,sys[i].params,NULL,0,sys[i].list->head->y_approx,globals->type,sys[i].diff_start,sys[i].no_ini_start,sys[i].user,NULL);
    //            custom_model->InitializeEqs(globals->global_params, sys[i].params, NULL, 0, sys[i].list->head->y_approx, globals->type, sys[i].diff_start, sys[i].no_ini_start, sys[i].user, NULL); //!!!! Should all states be reset? !!!!
    //    for (i = 0; i < asynch->globals->num_forcings; i++)
    //    {
    //        if (asynch->forcings[i]->flag == 3)	//!!!! I think .mon and binary files need this too !!!!
    //            Asynch_Set_Forcing_State(asynch, i, t_b, asynch->forcings[i]->first_file, asynch->forcings[i]->last_file);
    //    }
    //    for (i = 0; i < max_or_steps; i++)
    //    {
    //        globals->maxtime = t_b + (i)* inc;
    //        if (i)	Asynch_Advance(asynch, 0);

    //        for (j = 0; j < num_obs; j++)
    //        {
    //            owner = assignments[obs_locs[j]];
    //            my_link = (owner == my_rank);

    //            if (my_link)
    //                q[i*num_obs + j] = sys[obs_locs[j]]->list->tail->y_approx.ve[0];

    //            MPI_Bcast(&(q[i*num_obs + j]), 1, MPI_DOUBLE, owner, MPI_COMM_WORLD);
    //        }
    //    }

    //    if (my_rank == 0)
    //    {
    //        //printf("q after\n");
    //        //Print_VECTOR(q,num_total_obs);

    //        double second_diff = ComputeDiff(d_els, q, num_total_obs);
    //        printf("\nDifferences between q and data are %e %e\n\n", first_diff, second_diff);
    //    }
    //}



    //Clean up
    VecDestroy(&d);	//!!!! Blah !!!!

    MPI_Barrier(MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (asynch->my_rank == 0)
        printf("Total time for linear least squares fit: %.0f\n", stop - start);

    return 0;
    //if(second_diff < first_diff)	return 1;
    //else				return 0;
}


double ComputeDiff(const double * const d, const double * const q, unsigned int size)
{
    unsigned int i;
    double result = 0.0;

    for (i = 0; i < size; i++)
        result += (d[i] - q[i]) * (d[i] - q[i]);

    //return pow(result, 0.5);
    return result;
}


void ResetSysLS(Link* sys, unsigned int N, GlobalVars* globals, double t_0, double* x_start, unsigned int problem_dim, unsigned int num_forcings, TransData* my_data)
{
    unsigned i, j, k, l;
    Link* current;

    Flush_TransData(my_data);

    for (i = 0; i < N; i++)
    {
        current = &sys[i];
        if (current->my != NULL)
        {
            while (current->current_iterations > 1)
            {
                Remove_Head_Node(&current->my->list);
                (current->current_iterations)--;
            }
            current->my->list.head->t = t_0;
            current->last_t = t_0;
            current->steps_on_diff_proc = 1;
            current->iters_removed = 0;
            current->rejected = 0;
            if (current->num_parents == 0)
                current->ready = 1;
            else
                current->ready = 0;
            for (j = 0; j < problem_dim; j++)
                current->my->list.head->y_approx[j] = x_start[i*problem_dim + j];
            //v_copy(backup[i],current->my->list.head->y_approx);

            //Reset the next_save time
            if (current->save_flag)
            {
                current->next_save = t_0;		//!!!! This forces the print times to match up with the assimilation times !!!!
                                                //current->disk_iterations = 1;
            }

            //Reset peak flow information
            current->peak_time = t_0;
            dcopy(current->my->list.head->y_approx, current->peak_value, 0, current->dim);

            //Set hydrograph scale
            //current->Q_TM = backup[i]->ve[0];

            //Reset current state
            if (current->check_state != NULL)
                current->state = current->check_state(
                    current->my->list.head->y_approx, current->dim,
                    globals->global_params, globals->num_global_params,
                    current->params, globals->num_params,
                    current->qvs, current->has_dam, NULL);
            current->my->list.head->state = current->state;

            //Set forcings
            if (current->my->forcing_data)
            {
                for (k = 0; k < num_forcings; k++)
                {
                    if (current->my->forcing_values[k])
                    {
                        //Find the right index in forcings
                        for (l = 0; l < current->my->forcing_data[k].num_points - 1; l++)
                            if (current->my->forcing_data[k].data[l].time <= t_0 && t_0 < current->my->forcing_data[k].data[l + 1].time)	break;
                        double rainfall_buffer = current->my->forcing_data[k].data[l].value;
                        current->my->forcing_values[k] = rainfall_buffer;
                        current->my->forcing_indices[k] = l;

                        //Find and set the new change in data
                        for (j = l + 1; j < current->my->forcing_data[k].num_points; j++)
                        {
                            if (fabs(current->my->forcing_data[k].data[j].value - rainfall_buffer) > 1e-12)
                            {
                                current->my->forcing_change_times[k] = current->my->forcing_data[k].data[j].time;
                                break;
                            }
                        }
                        if (j == current->my->forcing_data[k].num_points)
                            current->my->forcing_change_times[k] = current->my->forcing_data[k].data[j - 1].time;

                        //Select new step size
                        //current->h = InitialStepSize(current->last_t,current,globals,workspace);
                    }
                }
            }
        }
    }
}

void Print_MATRIX(double** A, unsigned int m, unsigned int n)
{
    unsigned int i, j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)	printf("%.2e ", A[i][j]);
        printf(";\n");
    }
    printf("\n");
}

void Print_VECTOR(double* v, unsigned int dim)
{
    unsigned int i;

    for (i = 0; i < dim; i++)
        printf("[%d]: %.2e\n", i, v[i]);
    printf(";\n");
}

