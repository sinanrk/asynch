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

#include <stdbool.h>

#include "metis.h"

#include "assim_models.h"

int* Partition_METIS_ByEqs(Link* sys, unsigned int N, Link** leaves, unsigned int numleaves, unsigned int** my_sys, unsigned int* my_N, TransData* my_data, short int *getting)
{
    unsigned int i, j, start_index, end_index, loc, retval;
    unsigned int nodes_per_proc = numleaves / np;	//Number of leaves assigned to each process (except the last)
    Link* current;

    start_index = nodes_per_proc * my_rank;
    if (my_rank == np - 1)		end_index = numleaves;
    else				end_index = nodes_per_proc * (my_rank + 1);
    //	*my_N = end_index - start_index;
    *my_N = 0;
    unsigned int my_max_nodes = N - numleaves + nodes_per_proc;
    *my_sys = (unsigned int*)malloc(my_max_nodes * sizeof(unsigned int));	//The indices of this processes links (in sys)
    for (i = 0; i < my_max_nodes; i++)	(*my_sys)[i] = -1;
    for (i = start_index; i < end_index; i++)	(*my_sys)[i - start_index] = leaves[i]->location;
    for (i = 0; i < N; i++)	getting[i] = 0;

    //Initialize assignments
    int* assignments = (int*)malloc(N * sizeof(int));
    for (i = 0; i < N; i++)	assignments[i] = -1;

    //Form the graph to partition
    idx_t* xadj = malloc((N + 1) * sizeof(idx_t));
    idx_t* adjncy = malloc(2 * (N - 1) * sizeof(idx_t));
    idx_t index = 0;

    for (i = 0; i < N; i++)
    {
        xadj[i] = index;
        current = &sys[i];
        if (current->child != NULL)
        {
            adjncy[index] = current->child->location;
            index++;
        }
        for (j = 0; j < current->num_parents; j++)
        {
            adjncy[index] = current->parents[j]->location;
            index++;
        }
    }
    xadj[N] = 2 * (N - 1);

    //Partition the system
    idx_t nverts = N;
    idx_t parts = np;
    idx_t ncon = 1;
    idx_t objval;
    idx_t* partitions = (idx_t*)calloc(N, sizeof(idx_t));
    idx_t* vwgt = (idx_t*)malloc(N * sizeof(idx_t));

    for (i = 0; i < N; i++)
        vwgt[i] = ((UpstreamData*)(sys[i].user))->dim;
    if (np != 1)
    {
        retval = METIS_PartGraphKway(&nverts, &ncon, xadj, adjncy, vwgt, vwgt, NULL, &parts, NULL, NULL, NULL, &objval, partitions);
        //retval = METIS_PartGraphKway(&nverts,&ncon,xadj,adjncy,NULL,NULL,NULL,&parts,NULL,NULL,NULL,&objval,partitions);
        if (retval != METIS_OK)
        {
            printf("Error: METIS returned error code %i.\n", retval);
            return NULL;
        }
    }

    *my_N = 0;
    for (i = 0; i < N; i++)
    {
        assignments[i] = partitions[i];	//!!!! Just use assignments? !!!!
        if (partitions[i] == my_rank)
        {
            (*my_sys)[*my_N] = i;
            (*my_N)++;
        }
    }

    //Set the getting array and determine number of sending and receiving links
    for (i = 0; i < *my_N; i++)
    {
        //Receiving
        for (j = 0; j < sys[(*my_sys)[i]].num_parents; j++)
        {
            loc = sys[(*my_sys)[i]].parents[j]->location;
            if (assignments[loc] != my_rank)
            {
                getting[loc] = 1;
                my_data->receive_size[assignments[loc]]++;
            }
        }

        //Sending
        if (sys[(*my_sys)[i]].child != NULL)
        {
            loc = sys[(*my_sys)[i]].child->location;
            if (assignments[loc] != my_rank)
                my_data->send_size[assignments[loc]]++;
        }
    }

    //Reorder my_sys so that the links with lower numbering are towards the beginning
    merge_sort_distance(sys, *my_sys, *my_N);

    //Allocate space in my_data for recieving and sending
    for (j = 0; j < np; j++)
    {
        my_data->receive_data[j] = (Link**)malloc(my_data->receive_size[j] * sizeof(Link*));
        my_data->send_data[j] = (Link**)malloc(my_data->send_size[j] * sizeof(Link*));
    }

    //Set the receive_data and send_data arrays
    int* current_receive_size = (int*)calloc(np, sizeof(int));
    int* current_send_size = (int*)calloc(np, sizeof(int));
    for (i = 0; i < *my_N; i++)
    {
        //Receiving
        for (j = 0; j < sys[(*my_sys)[i]].num_parents; j++)
        {
            loc = sys[(*my_sys)[i]].parents[j]->location;
            if (assignments[loc] != my_rank)
            {
                my_data->receive_data[assignments[loc]][current_receive_size[assignments[loc]]] = &sys[loc];
                current_receive_size[assignments[loc]]++;
            }
        }

        //Sending
        if (sys[(*my_sys)[i]].child != NULL)
        {
            loc = sys[(*my_sys)[i]].child->location;
            if (assignments[loc] != my_rank)
            {
                my_data->send_data[assignments[loc]][current_send_size[assignments[loc]]] = &sys[(*my_sys)[i]];
                current_send_size[assignments[loc]]++;
            }
        }
    }

    //Clean up
    free(current_receive_size);
    free(current_send_size);
    free(xadj);
    free(adjncy);
    free(partitions);
    free(vwgt);
    (*my_sys) = (unsigned int*)realloc(*my_sys, *my_N * sizeof(unsigned int));

    return assignments;
}


void Setup_Errors(AsynchSolver* asynch, unsigned int problem_dim)
{
    GlobalVars* GlobalVars = asynch->GlobalVars;
    ErrorData* GlobalErrors = asynch->GlobalErrors;
    unsigned int i, max_dim = GlobalVars->max_dim;

    GlobalErrors->abstol.ve = realloc(GlobalErrors->abstol.ve, max_dim * sizeof(double));
    GlobalErrors->reltol.ve = realloc(GlobalErrors->reltol.ve, max_dim * sizeof(double));
    GlobalErrors->abstol_dense.ve = realloc(GlobalErrors->abstol_dense.ve, max_dim * sizeof(double));
    GlobalErrors->reltol_dense.ve = realloc(GlobalErrors->reltol_dense.ve, max_dim * sizeof(double));
    GlobalErrors->abstol.dim = GlobalErrors->reltol.dim = GlobalErrors->reltol_dense.dim = GlobalErrors->reltol_dense.dim = max_dim;

    //Setup error
    for (i = problem_dim + 1; i < max_dim; i++)
    {
        GlobalErrors->abstol.ve[i] = GlobalErrors->abstol.ve[problem_dim];
        GlobalErrors->reltol.ve[i] = GlobalErrors->reltol.ve[problem_dim];
        GlobalErrors->abstol_dense.ve[i] = GlobalErrors->abstol_dense.ve[problem_dim];
        GlobalErrors->reltol_dense.ve[i] = GlobalErrors->reltol_dense.ve[problem_dim];
    }
}

//Creates an array vareq_shift with allstates = sum(assim_dim).
//This is used to remove the sensitivities that are not used.
unsigned int BuildStateShift(AsynchSolver* asynch, unsigned int allstates, unsigned int* obs_locs, unsigned int num_obs, unsigned int** vareq_shift, unsigned int** inv_vareq_shift)
{
    unsigned int i, j, allstates_needed;
    Link *sys = asynch->sys;
    int *assignments = asynch->assignments;
    bool *needed = (bool*)calloc(allstates, sizeof(bool));

    for (i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == my_rank)
        {
            UpstreamData * updata = (UpstreamData*)sys[obs_locs[i]].user;
            for (j = 0; j < updata->num_fit_states; j++)
                needed[updata->fit_to_universal[j]] = true;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, needed, allstates, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

    *vareq_shift = (unsigned int*)malloc(allstates * sizeof(unsigned int));
    *inv_vareq_shift = (unsigned int*)calloc(allstates, sizeof(unsigned int));

    //Compute the shifts
    unsigned int shift = 0;
    for (i = 0; i < allstates; i++)
    {
        if (needed[i] == false)
            shift++;
        else
        {
            (*vareq_shift)[i] = i - shift;
            (*inv_vareq_shift)[i - shift] = i;
        }
    }

    /*
    printf("Shift is (%u, %u)\n",allstates,allstates-shift);
    for(i=0;i<allstates;i++)
    printf("%u ",(*vareq_shift)[i]);
    printf("\n");
    printf("Inv shift is (%u)\n",allstates-shift);
    for(i=0;i<allstates;i++)
    printf("%u ",(*inv_vareq_shift)[i]);
    printf("\n");
    */
    allstates_needed = allstates - shift;
    *inv_vareq_shift = (unsigned int*)realloc(*inv_vareq_shift, allstates_needed * sizeof(unsigned int));
    free(needed);
    return allstates_needed;
}

//Data assimilation model (Old Model 315) ************************************************************************************


void SetParamSizes_Assim(GlobalVars* GlobalVars, void* external)
{
    GlobalVars->uses_dam = 0;
    GlobalVars->params_size = 20;
    GlobalVars->dam_params_size = 0;
    GlobalVars->area_idx = 2;
    GlobalVars->areah_idx = 1;
    GlobalVars->disk_params = 12;
    GlobalVars->convertarea_flag = 0;
    GlobalVars->num_forcings = 1;
    GlobalVars->min_error_tolerances = 4;
}


void ConvertParams_Assim(VEC params, unsigned int type, void* external)
{
    params.ve[0] *= 1000;	//L: km -> m
    params.ve[3] *= .001;	//h_b: mm -> m
    params.ve[4] *= .001;	//h_H: mm -> m
}

void InitRoutines_Assim(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    UpstreamData* updata = (UpstreamData*)(link->user);
    unsigned int i, problem_dim = 2;	//Number of model eqs

    link->dim = problem_dim + problem_dim + (problem_dim - 1)*(problem_dim - 1) //Model eqs + variational eqs from this link
        + updata->num_upstreams * problem_dim;	                                //Variational eqs from upstreams !!!! Too high? !!!!
    //for(i=0;i<link->num_parents;i++)
    //	link->dim += updata->num_upstreams[i] * problem_dim;	//Variational eqs from upstreams !!!! Too high? !!!!

    if (link->dim != updata->dim)
        printf("[%i]: Warning: calculated the number of equations at link %u twice and got %u and %u.\n", my_rank, link->ID, link->dim, updata->dim);
    link->no_ini_start = 2;
    link->diff_start = 0;

    link->num_dense = link->dim - 1;	//Take out s_p
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;
    for (i = 1; i < link->num_dense; i++)	link->dense_indices[i] = i + 1;

    link->f = &assim_river_rainfall_adjusted_custom;
    link->alg = NULL;
    link->state_check = NULL;
    link->CheckConsistency = &CheckConsistency_Nonzero_2States;
    link->RKSolver = &ExplicitRKSolver;
}


void InitRoutines_Model(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    UpstreamData* updata = (UpstreamData*)(link->user);
    unsigned int problem_dim = 2;	//Number of model eqs

    link->dim = 2;
    link->no_ini_start = 2;
    link->diff_start = 0;

    link->num_dense = 1;
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;

    link->f = &river_rainfall_adjusted;
    link->alg = NULL;
    link->state_check = NULL;
    link->CheckConsistency = &CheckConsistency_Nonzero_2States;
    link->RKSolver = &ExplicitRKSolver;
}

void Precalculations_Assim(Link* link_i, VEC global_params, VEC params, unsigned int disk_params, unsigned int params_size, unsigned short int dam, unsigned int type, void* external)
{
    //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
    //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
    //Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
    //The numbering is:        0      1        2     3   4   5
    //Need to set entries 12-19 of params.
    double* vals = params.ve;
    double K_T = 1.0;
    double s_r = 1.0;
    double rootS_h = pow(vals[7], .5);
    double L = params.ve[0];
    double A_h = params.ve[1] * 1e6;	//Put into m^2
    double eta = params.ve[8];
    double v_r = global_params.ve[0];
    double lambda_1 = global_params.ve[1];
    double lambda_2 = global_params.ve[2];
    double v_h = global_params.ve[3];
    double A_r = global_params.ve[4];
    double RC = global_params.ve[5];

    //!!!! Clean this model. You don't really need 20 parameters... !!!!
    vals[12] = 60.0*v_r*pow(vals[2] / A_r, lambda_2) / ((1.0 - lambda_1)*vals[0]);	//invtau [1/min]
    vals[13] = vals[3] / s_r; //epsilon
    vals[14] = v_h*L;	//c_1 [m^2/s]
    vals[15] = vals[6] * vals[0] * vals[3] / 3600.0; //c_2
    vals[16] = (1e-3 / 60.0) * RC;	//c_3
    vals[17] = 60.0*v_h*L / A_h;	//c_4 [1/min], A_h converted above
    vals[18] = K_T / 60.0;
    vals[19] = vals[6] / (60.0*s_r);

    //iparams.ve[0] = link_i->location; //!!!! Is this even needed anywhere? !!!!
}

int ReadInitData_Assim(VEC global_params, VEC params, QVSData* qvs, unsigned short int dam, VEC y_0, unsigned int type, unsigned int diff_start, unsigned int no_init_start, void* user, void* external)
{
    //For this type, all initial conditions for variational equation must be set here.
    //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
    //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
    //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
    //The numbering is:        0      1        2     3   4   5
    unsigned int i;
    unsigned int offset = 2;

    y_0.ve[offset] = 1.0;
    y_0.ve[offset + 1] = 1.0;
    y_0.ve[offset + 2] = 0.0;
    for (i = offset + 3; i < y_0.dim; i++)	y_0.ve[i] = 0.0;

    return 0;
}

//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void assim_river_rainfall_adjusted_custom(double t, VEC y_i, VEC* y_p, unsigned short int num_parents, VEC global_params, double* forcing_values, QVSData* qvs, VEC params, int state, void* user, VEC ans)
{
    unsigned int i, j;
    unsigned int dim = ans.dim;
    unsigned int offset = 2;		//!!!! This needs to be num_dense, but without variational eqs !!!!
    unsigned int parent_offset;
    unsigned int problem_dim = 2;
    unsigned int all_states = (dim - offset) / problem_dim;
    double inflow = 0.0;
    UpstreamData* updata = (UpstreamData*)user;

    double q = y_i.ve[0];
    double s_p = y_i.ve[1];

    double L = params.ve[0];
    double invtau = params.ve[12];
    double c_1 = params.ve[14];
    double c_3 = params.ve[16];
    double c_4 = params.ve[17];
    double lambda_1 = global_params.ve[1];

    double q_to_lambda_1 = pow(q, lambda_1);
    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);
    double deriv_qpl = 1.0;

    double q_pl = s_p;

    //Flux equation (y_i[0])
    ans.ve[0] = -q + c_1 * q_pl;
    for (i = 0; i < num_parents; i++)
        inflow += y_p[i].ve[0];
    ans.ve[0] = invtau * q_to_lambda_1 * (inflow + ans.ve[0]);

    //Ponded water equation (y_i[1])
    ans.ve[1] = c_3 * forcing_values[0] - c_4 * q_pl;
    //ans.ve[1] = c_3 * ( max(forcing_values[0] + 20.0*sin(t/5.0),0.0)) - c_4 * q_pl;

    //!!!! Pull if statements out of loops (should just need two cases total) !!!!
    //!!!! A lot of terms get repeated !!!!

    //Eqs for variational equations
    for (i = offset; i < dim; i++)	ans.ve[i] = 0.0;

    //s variable from this link
    ans.ve[offset] = -c_4*deriv_qpl*y_i.ve[offset];

    //q variables from this link
//	if(lambda_1 > 1e-12 && (inflow) > 1e-12)
    ans.ve[offset + 1] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i.ve[offset + 1];
    //	else
    //		ans.ve[offset + 1] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i.ve[offset + 1];

    //	if(lambda_1 > 1e-12 && (inflow) > 1e-12)
    ans.ve[offset + 2] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i.ve[offset + 2] + invtau*c_1*q_to_lambda_1*deriv_qpl * y_i.ve[offset];
    //	else
    //		ans.ve[offset + 2] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i.ve[offset + 2] + invtau*c_1*deriv_qpl*y_i.ve[offset];

        //Adjust offset
    offset += 3;

    //Variables from parents
    for (i = 0; i < num_parents; i++)
    {
        parent_offset = 1 + problem_dim;

        //Get the number of upstreams link for that parent
        unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;

        for (j = 0; j < num_upstreams; j++)
        {
            ans.ve[offset] = invtau * q_to_lambda_1 * y_p[i].ve[parent_offset];
            //			if(lambda_1 > 1e-12 && (inflow) > 1e-12)
            ans.ve[offset] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i.ve[offset];
            //			else
            //				ans.ve[offset] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i.ve[offset];

            ans.ve[offset + 1] = invtau * q_to_lambda_1 * y_p[i].ve[parent_offset + 1];
            //			if(lambda_1 > 1e-12 && (inflow) > 1e-12)
            ans.ve[offset + 1] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i.ve[offset + 1];
            //			else
            //				ans.ve[offset + 1] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i.ve[offset + 1];

            offset += 2;
            parent_offset += 2;
        }
    }
}



//Data assimilation model (Model 254) ************************************************************************************


void SetParamSizes_Assim_254(GlobalVars* GlobalVars, void* external)
{
    GlobalVars->uses_dam = 0;
    GlobalVars->params_size = 8;
    GlobalVars->dam_params_size = 0;
    GlobalVars->area_idx = 0;
    GlobalVars->areah_idx = 2;
    GlobalVars->disk_params = 3;
    GlobalVars->convertarea_flag = 0;
    GlobalVars->num_forcings = 3;
    GlobalVars->min_error_tolerances = 8;
}


void ConvertParams_Assim_254(VEC params, unsigned int type, void* external)
{
    params.ve[1] *= 1000;		//L_h: km -> m
    params.ve[2] *= 1e6;		//A_h: km^2 -> m^2
}

void InitRoutines_Assim_254(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    UpstreamData* updata = (UpstreamData*)(link->user);
    unsigned int i, problem_dim = 4;	//Number of model eqs

    link->dim = (problem_dim + 3) + 11	    //Model eqs + variational eqs from this link
        + updata->num_upstreams * problem_dim;	//Variational eqs from upstreams
    //for(i=0;i<link->num_parents;i++)
    //	link->dim += updata->num_upstreams[i] * problem_dim;	//Variational eqs from upstreams
    //if(link->dim != updata->dim)
    //	printf("[%i]: Warning: calculated the number of equations at link %u twice and got %u and %u.\n",my_rank,link->ID,link->dim,updata->dim);
    link->no_ini_start = 4;
    link->diff_start = 0;

    link->num_dense = link->dim - 5;	//Only q, q_b, variational eqs
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;
    link->dense_indices[1] = 6;
    for (i = 2; i < link->num_dense; i++)	link->dense_indices[i] = i + 5;

    link->f = &TopLayerHillslope_extras_assim;
    link->alg = NULL;
    link->state_check = NULL;
    link->CheckConsistency = &CheckConsistency_Nonzero_Model254;
    link->RKSolver = &ExplicitRKSolver;
}

void InitRoutines_Model_254(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    link->dim = 7;
    link->no_ini_start = 4;
    link->diff_start = 0;

    link->num_dense = 2;
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;
    link->dense_indices[1] = 6;

    if (link->res)
    {
        link->f = &TopLayerHillslope_Reservoirs;
        link->RKSolver = &ForcedSolutionSolver;
    }
    else			link->f = &TopLayerHillslope_extras;
    link->alg = NULL;
    link->state_check = NULL;
    link->CheckConsistency = &CheckConsistency_Nonzero_AllStates_q;
}

void Precalculations_Assim_254(Link* link_i, VEC global_params, VEC params, unsigned int disk_params, unsigned int params_size, unsigned short int dam, unsigned int type, void* external)
{
    //Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
    //The numbering is:	0   1   2    3     4   5   6   7 
    //Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
    //The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
    double* vals = params.ve;
    double A_i = params.ve[0];
    double L_i = params.ve[1];
    double A_h = params.ve[2];

    double v_0 = global_params.ve[0];
    double lambda_1 = global_params.ve[1];
    double lambda_2 = global_params.ve[2];
    double v_h = global_params.ve[3];
    double k_i_factor = global_params.ve[5];

    vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
    vals[4] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
    vals[5] = vals[4] * k_i_factor;	//[1/min] k_i
    vals[6] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
    vals[7] = A_h / 60.0;	//  c_2
}

int ReadInitData_Assim_254(VEC global_params, VEC params, QVSData* qvs, unsigned short int dam, VEC y_0, unsigned int type, unsigned int diff_start, unsigned int no_init_start, void* user, void* external)
{
    //For this type, all initial conditions for variational equation must be set here.
    unsigned int i;
    unsigned int offset = 7;

    //For this type, the extra states need to be set (4,5,6)
    y_0.ve[4] = 0.0;
    y_0.ve[5] = 0.0;
    y_0.ve[6] = y_0.ve[0];

    y_0.ve[offset++] = 1.0;  //ds_p/ds_p0
    y_0.ve[offset++] = 0.0;  //ds_p/ds_t0

    y_0.ve[offset++] = 0.0;  //ds_t/ds_p0
    y_0.ve[offset++] = 1.0;  //ds_t/ds_t0

    y_0.ve[offset++] = 0.0;  //ds_s/ds_p0
    y_0.ve[offset++] = 0.0;  //ds_s/ds_t0
    y_0.ve[offset++] = 1.0;  //ds_s/ds_s0

    y_0.ve[offset++] = 1.0;  //dq/dq_0
    y_0.ve[offset++] = 0.0;  //dq/ds_p0
    y_0.ve[offset++] = 0.0;  //dq/ds_t0
    y_0.ve[offset++] = 0.0;  //dq/ds_s0

    for (i = offset; i < y_0.dim; i++)	y_0.ve[i] = 0.0;	//From upstreams

    return 0;
}

void CheckConsistency_Nonzero_Model254(VEC y, VEC params, VEC global_params)
{
    unsigned int i, problem_dim = 7;

    if (y.ve[0] < 1e-14)	y.ve[0] = 1e-14;
    for (i = 1; i < problem_dim; i++)
        if (y.ve[i] < 0.0)	y.ve[i] = 0.0;
}

//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: A_i,L_i,A_h | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2       3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void TopLayerHillslope_extras_assim(double t, VEC y_i, VEC* y_p, unsigned short int num_parents, VEC global_params, double* forcing_values, QVSData* qvs, VEC params, int state, void* user, VEC ans)
{
    unsigned int i, j;

    UpstreamData *updata = (UpstreamData*)user;

    double lambda_1 = global_params.ve[1];
    double k_3 = global_params.ve[4];	//[1/min]
    double h_b = global_params.ve[6];	//[m]
    double S_L = global_params.ve[7];	//[m]
    double A = global_params.ve[8];
    double B = global_params.ve[9];
    double exponent = global_params.ve[10];
    double v_B = global_params.ve[11];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params.ve[1];	//[m]
    double A_h = params.ve[2];	//[m^2]
    double invtau = params.ve[3];	//[1/min]
    double k_2 = params.ve[4];	//[1/min]
    double k_i = params.ve[5];	//[1/min]
    double c_1 = params.ve[6];
    double c_2 = params.ve[7];

    double q = y_i.ve[0];		//[m^3/s]
    double s_p = y_i.ve[1];	//[m]
    double s_t = y_i.ve[2];	//[m]
    double s_s = y_i.ve[3];	//[m]
    //double s_precip = y_i.ve[4];	//[m]
    //double V_r = y_i.ve[5];	//[m^3]
    double q_b = y_i.ve[6];	//[m^3/s]

    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    //A few calculations...
    double q_to_lambda_1 = pow(q, lambda_1);
    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);
    double remaining = 1.0 - s_t / S_L;
    double pow_term = (remaining > 0.0) ? pow(remaining, exponent) : 0.0;
    double pow_term_m1 = (remaining > 1e-12) ? pow_term / remaining : pow(1e-12, exponent - 1.0);
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

    //Discharge
    double inflow = 0.0;
    ans.ve[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i < num_parents; i++)
        inflow += y_p[i].ve[0];
    ans.ve[0] = invtau * q_to_lambda_1 * (inflow + ans.ve[0]);

    //Hillslope
    ans.ve[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans.ve[2] = q_pt - q_ts - e_t;
    ans.ve[3] = q_ts - q_sl - e_s;

    //Additional states
    ans.ve[4] = forcing_values[0] * c_1;
    ans.ve[5] = q_pl;
    ans.ve[6] = q_sl * A_h - q_b*60.0;
    for (i = 0; i < num_parents; i++)
        ans.ve[6] += y_p[i].ve[6] * 60.0;
    //ans.ve[6] += k_3*y_p[i].ve[3]*A_h;
    ans.ve[6] *= v_B / L;


    //!!!! Pull if statements out of loops (should just need two cases total) !!!!
    //!!!! A lot of terms get repeated !!!!

    //Init for variational equations
    unsigned int offset = 7, dim = ans.dim, problem_dim = 7;
    for (i = offset; i < dim; i++)	ans.ve[i] = 0.0;

    //Compute partial derivatives (local variables)
    double dfq_dq = lambda_1 * invtau * q_to_lambda_1_m1 * (-q + c_2*(k_2*s_p + k_3*s_s) + inflow) - invtau * q_to_lambda_1;
    double dfq_dsp = invtau*q_to_lambda_1*c_2*k_2;
    double dfq_dss = invtau*q_to_lambda_1*c_2*k_3;
    double dfsp_dsp = -k_2 - k_t;
    double dfsp_dst = k_2 / S_L*B*exponent*pow_term_m1*s_p;
    double dfst_dsp = k_t;
    double dfst_dst = -dfsp_dst - k_i;
    double dfss_dst = k_i;
    double dfss_dss = -k_3;

    //Compute partial derivatives (upstreams variables)
    double dfq_dupq = invtau*q_to_lambda_1;

    //Hillslope variational eqs
    //!!!! I think these are needed only if changing the appropriate state... !!!!
    ans.ve[offset] = dfsp_dsp * y_i.ve[offset] + dfsp_dst * y_i.ve[offset + 2];	//s_p, s_p
    ans.ve[offset + 1] = dfsp_dsp * y_i.ve[offset + 1] + dfsp_dst * y_i.ve[offset + 3];	//s_p, s_t
    ans.ve[offset + 2] = dfst_dsp * y_i.ve[offset] + dfst_dst * y_i.ve[offset + 2];	//s_t, s_p
    ans.ve[offset + 3] = dfst_dsp * y_i.ve[offset + 1] + dfst_dst * y_i.ve[offset + 3];	//s_t, s_t	
    ans.ve[offset + 4] = dfss_dst * y_i.ve[offset + 2] + dfss_dss * y_i.ve[offset + 4];	//s_s, s_p
    ans.ve[offset + 5] = dfss_dst * y_i.ve[offset + 3] + dfss_dss * y_i.ve[offset + 5];	//s_s, s_t
    ans.ve[offset + 6] = dfss_dss * y_i.ve[offset + 6];	//s_s, s_s

    //Discharge variational eqs from this link
    ans.ve[offset + 7] = dfq_dq * y_i.ve[offset + 7]; //q, q
    ans.ve[offset + 8] = dfq_dq * y_i.ve[offset + 8] + dfq_dsp * y_i.ve[offset + 0] + dfq_dss * y_i.ve[offset + 4]; //q, s_p
    ans.ve[offset + 9] = dfq_dq * y_i.ve[offset + 9] + dfq_dsp * y_i.ve[offset + 1] + dfq_dss * y_i.ve[offset + 5]; //q, s_t
    ans.ve[offset + 10] = dfq_dq * y_i.ve[offset + 10] + dfq_dss * y_i.ve[offset + 6]; //q, s_s

    //Discharge variational eqs from parent links
    unsigned int current_idx = offset + 11, parent_idx;
    for (i = 0; i < num_parents; i++)
    {
        //parent_idx = offset + 8;
        parent_idx = offset + 7;

        //Get the number of upstreams link for that parent
        unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;

        for (j = 0; j < num_upstreams; j++)
        {
            ans.ve[current_idx] = dfq_dupq * y_p[i].ve[parent_idx] + dfq_dq * y_i.ve[current_idx]; //q, upq
            ans.ve[current_idx + 1] = dfq_dupq * y_p[i].ve[parent_idx + 1] + dfq_dq * y_i.ve[current_idx + 1]; //q, ups_p
            ans.ve[current_idx + 2] = dfq_dupq * y_p[i].ve[parent_idx + 2] + dfq_dq * y_i.ve[current_idx + 2]; //q, ups_t
            ans.ve[current_idx + 3] = dfq_dupq * y_p[i].ve[parent_idx + 3] + dfq_dq * y_i.ve[current_idx + 3]; //q, ups_s
            current_idx += 4;
            parent_idx += 4;
        }
    }
}

//I think these values need to be set only for assigned links, not getting.
//!!!! This data is really only needed at the gauges. !!!!
//fit_states[i] holds the index in each state vector of the ith sensitivity at this link.
//fit_to_universal[i] holds universal index of the ith sensitivity at this link.
//These only store those sensitivites used for the fitting.
void Setup_Fitting_Data_Model254(AsynchSolver* asynch, unsigned int* obs_locs, unsigned int num_obs)
{
    unsigned int i, j, k, my_N = asynch->my_N, *my_sys = asynch->my_sys, *assignments = asynch->assignments;
    Link *sys = asynch->sys, *current;
    UpstreamData *updata;

    //Number of states to fit
    unsigned int counter;
    unsigned int problem_dim = 7;	//!!!! Should be allowed to vary by link !!!!
    unsigned int assim_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    //unsigned int num_change_states = 1;	//For q
    unsigned int num_change_states = 2;	//For q and s_p

    //Find links upstreams from gauges
    bool *is_above_gauges;
    unsigned int *above_gauges;
    unsigned int num_above = GaugeDownstream(asynch, obs_locs, num_obs, &above_gauges, &is_above_gauges);

    //Calculate the number of states needed for the fitting
    //unsigned int allstates_needed = num_above;	//For q
    //unsigned int allstates_needed = num_above * 2;	//For q and s_p

    //for(i=0;i<my_N;i++)
    for (i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == my_rank)
        {
            //current = sys[my_sys[i]];
            current = &sys[obs_locs[i]];
            updata = (UpstreamData*)current->user;

            updata->num_fit_states = 1;
            for (j = 0; j < current->num_parents; j++)
            {
                //Get the number of upstreams link for that parent
                unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
                updata->num_fit_states += num_upstreams;
            }
            updata->num_fit_states *= num_change_states;

            updata->fit_states = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));
            updata->fit_to_universal = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));
            /*
                        //For q
                        updata->fit_states[0] = problem_dim + 7;	//q, q here
                        updata->fit_to_universal[0] = current->location * assim_dim;
                        counter = 1;
                        for(j=0;j<current->num_parents;j++)
                        {
                            for(k=0;k<updata->num_upstreams[j];k++)
                            {
                                updata->fit_states[counter] = problem_dim + 7 + assim_dim*counter;
                                updata->fit_to_universal[counter] = updata->upstreams[j][k] * assim_dim;
                                counter++;
                            }
                        }
            */

            //For q, s_p
            updata->fit_states[0] = problem_dim + 7;	//q, q here
            updata->fit_to_universal[0] = current->location * assim_dim;
            updata->fit_states[1] = problem_dim + 8;	//q, s_p here
            updata->fit_to_universal[1] = current->location * assim_dim + 1;
            counter = 2;
            for (j = 0; j < current->num_parents; j++)
            {
                unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
                Link **upstreams = ((UpstreamData *)updata->upstreams[i]->user)->upstreams;

                for (k = 0; k < num_upstreams; k++)
                {
                    updata->fit_states[counter] = problem_dim + 7 + assim_dim*counter / 2;
                    updata->fit_to_universal[counter] = upstreams[k]->location * assim_dim;
                    updata->fit_states[counter + 1] = problem_dim + 8 + assim_dim*counter / 2;
                    updata->fit_to_universal[counter + 1] = upstreams[k]->location * assim_dim + 1;
                    counter += 2;
                }
            }

        }
    }


    free(is_above_gauges);
    free(above_gauges);
    //return allstates_needed;
}


//Data assimilation model (Model 254 trimmed, q) ************************************************************************************


//For modifying q
void InitRoutines_Assim_254_q(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    UpstreamData* updata = (UpstreamData*)(link->user);
    unsigned int i, problem_dim = 4;	//Number of model eqs

    //For q only
    link->dim = problem_dim + 1	//Model eqs + variational eqs from this link
        + updata->num_upstreams; //Variational eqs from upstreams
    //for(i=0;i<link->num_parents;i++)
    //	link->dim += updata->num_upstreams[i];	//Variational eqs from upstreams
    link->no_ini_start = 4;
    link->diff_start = 0;

    link->num_dense = link->dim - 3;	//Only q, variational eqs
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;
    for (i = 1; i < link->num_dense; i++)	link->dense_indices[i] = i + 3;

    link->f = &TopLayerHillslope_assim_q;
    link->alg = NULL;
    link->state_check = NULL;
    link->CheckConsistency = &CheckConsistency_Nonzero_Model252;
    link->RKSolver = &ExplicitRKSolver;
}

void InitRoutines_Model_252(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    link->dim = 4;
    link->no_ini_start = 4;
    link->diff_start = 0;

    link->num_dense = 1;
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;

    if (link->res)
    {
        link->f = &TopLayerHillslope_Reservoirs;
        link->RKSolver = &ForcedSolutionSolver;
    }
    else			link->f = &TopLayerHillslope;
    link->alg = NULL;
    link->state_check = NULL;
    link->CheckConsistency = &CheckConsistency_Nonzero_AllStates_q;
}


int ReadInitData_Assim_254_q(VEC global_params, VEC params, QVSData* qvs, unsigned short int dam, VEC y_0, unsigned int type, unsigned int diff_start, unsigned int no_init_start, void* user, void* external)
{
    //For this type, all initial conditions for variational equation must be set here.
    unsigned int i;
    unsigned int offset = 4;

    y_0.ve[offset++] = 1.0;  //dq/dq_0

    for (i = offset; i < y_0.dim; i++)
        y_0.ve[i] = 0.0;	//From upstreams

    return 0;
}

void CheckConsistency_Nonzero_Model252(VEC y, VEC params, VEC global_params)
{
    unsigned int i, problem_dim = 4;

    if (y.ve[0] < 1e-14)	y.ve[0] = 1e-14;
    for (i = 1; i < problem_dim; i++)
        if (y.ve[i] < 0.0)	y.ve[i] = 0.0;
}

//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: A_i,L_i,A_h | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2       3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10     
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void TopLayerHillslope_assim_q(double t, VEC y_i, VEC* y_p, unsigned short int num_parents, VEC global_params, double* forcing_values, QVSData* qvs, VEC params, int state, void* user, VEC ans)
{
    unsigned int i;

    UpstreamData *updata = (UpstreamData*)user;

    double lambda_1 = global_params.ve[1];
    double k_3 = global_params.ve[4];	//[1/min]
    double h_b = global_params.ve[6];	//[m]
    double S_L = global_params.ve[7];	//[m]
    double A = global_params.ve[8];
    double B = global_params.ve[9];
    double exponent = global_params.ve[10];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params.ve[1];	//[m]
    double A_h = params.ve[2];	//[m^2]
    double invtau = params.ve[3];	//[1/min]
    double k_2 = params.ve[4];	//[1/min]
    double k_i = params.ve[5];	//[1/min]
    double c_1 = params.ve[6];
    double c_2 = params.ve[7];

    double q = y_i.ve[0];		//[m^3/s]
    double s_p = y_i.ve[1];	//[m]
    double s_t = y_i.ve[2];	//[m]
    double s_s = y_i.ve[3];	//[m]

    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    //A few calculations...
    double q_to_lambda_1 = pow(q, lambda_1);
    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);
    double remaining = 1.0 - s_t / S_L;
    double pow_term = (remaining > 0.0) ? pow(remaining, exponent) : 0.0;
    double pow_term_m1 = (remaining > 1e-12) ? pow_term / remaining : pow(1e-12, exponent - 1.0);
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

    //Discharge
    double inflow = 0.0;
    ans.ve[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i < num_parents; i++)
        inflow += y_p[i].ve[0];
    ans.ve[0] = invtau * q_to_lambda_1 * (inflow + ans.ve[0]);

    //Hillslope
    ans.ve[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans.ve[2] = q_pt - q_ts - e_t;
    ans.ve[3] = q_ts - q_sl - e_s;


    //Init for variational equations
    unsigned int offset = 4, dim = ans.dim, problem_dim = 4;
    for (i = offset; i < dim; i++)	ans.ve[i] = 0.0;	//!!!! Is this needed? !!!!

    //Compute partial derivatives (local variables)
    double dfq_dq = lambda_1 * invtau * q_to_lambda_1_m1 * (-q + c_2*(k_2*s_p + k_3*s_s) + inflow) - invtau * q_to_lambda_1;
    double dfq_dsp = invtau*q_to_lambda_1*c_2*k_2;
    double dfq_dss = invtau*q_to_lambda_1*c_2*k_3;
    double dfsp_dsp = -k_2 - k_t;
    double dfsp_dst = k_2 / S_L*B*exponent*pow_term_m1*s_p;
    double dfst_dsp = k_t;
    double dfst_dst = -dfsp_dst - k_i;
    double dfss_dst = k_i;
    double dfss_dss = -k_3;

    //Compute partial derivatives (upstreams variables)
    double dfq_dupq = invtau*q_to_lambda_1;

    //Discharge variational eqs from this link
    ans.ve[offset] = dfq_dq * y_i.ve[offset]; //q, q

    //Discharge variational eqs from parent links
    unsigned int current_idx = offset + 1;

    for (i = 0; i < num_parents; i++)
    {
        unsigned int num_upstreams = ((UpstreamData *)updata->parents[i]->user)->num_upstreams;
        for (unsigned int j = 0; j < num_upstreams; j++)
        {
            unsigned int parent_idx = j + offset;

            ans.ve[current_idx] = dfq_dupq * y_p[i].ve[parent_idx] + dfq_dq * y_i.ve[current_idx]; //q, upq
            current_idx += 1;
        }
    }
}

//I think these values need to be set only for assigned links, not getting.
//!!!! This data is really only needed at the gauges. !!!!
//fit_states[i] holds the index in each state vector of the ith sensitivity at this link.
//fit_to_universal[i] holds universal index of the ith sensitivity at this link.
//These only store those sensitivites used for the fitting.
void Setup_Fitting_Data_Model254_q(AsynchSolver* asynch, unsigned int* obs_locs, unsigned int num_obs)
{
    unsigned int i, j, k, my_N = asynch->my_N, *my_sys = asynch->my_sys, *assignments = asynch->assignments;
    Link *sys = asynch->sys, *current;
    UpstreamData *updata;

    //Number of states to fit
    unsigned int counter;
    unsigned int problem_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int assim_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int num_change_states = 1;	//For q
    //unsigned int num_change_states = 2;	//For q and s_p

    ////Find links upstreams from gauges
    //short int* is_above_gauges;
    //unsigned int* above_gauges;
    //unsigned int num_above = GaugeDownstream(asynch,&above_gauges,&is_above_gauges,obs_locs,num_obs);

    //Calculate the number of states needed for the fitting
    //unsigned int allstates_needed = num_above * 2;	//For q and s_p
    //unsigned int allstates_needed = num_above;	//For q

    for (i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == my_rank)
        {
            current = &sys[obs_locs[i]];
            updata = (UpstreamData*)current->user;

            updata->num_fit_states = 1 + updata->num_upstreams;
            //for(j=0;j<current->num_parents;j++)
            //	updata->num_fit_states += updata->num_upstreams[j];
            updata->num_fit_states *= num_change_states;

            updata->fit_states = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));
            updata->fit_to_universal = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));

            //For q
            updata->fit_states[0] = problem_dim;	//q, q here
            updata->fit_to_universal[0] = current->location * assim_dim;
            counter = 1;
            //for(j=0;j<current->num_parents;j++)
            //{
            unsigned int num_upstreams = ((UpstreamData *)updata)->num_upstreams;
            Link **upstreams = ((UpstreamData *)updata)->upstreams;

            for (k = 0; k < num_upstreams; k++)
            {
                updata->fit_states[counter] = problem_dim + counter;
                updata->fit_to_universal[counter] = upstreams[k]->location * assim_dim;
                counter++;
            }
            //}


        }
    }

    //free(is_above_gauges);
    //free(above_gauges);
    //return allstates_needed;
}




//Data assimilation model (Model 254 trimmed, q and s_p) ************************************************************************************


//For modifying q and s_p
void InitRoutines_Assim_254_qsp(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    UpstreamData* updata = (UpstreamData*)(link->user);
    unsigned int i, problem_dim = 4;	//Number of model eqs

    //For q only
    link->dim = problem_dim + 5         //Model eqs + variational eqs from this link
        + updata->num_upstreams * 2;	    //Variational eqs from upstreams
    //for(i=0;i<link->num_parents;i++)
    //	link->dim += updata->num_upstreams[i] * 2;	
    link->no_ini_start = 4;
    link->diff_start = 0;

    link->num_dense = link->dim - 3;	//Only q, variational eqs
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;
    for (i = 1; i < link->num_dense; i++)	link->dense_indices[i] = i + 3;

    link->f = &TopLayerHillslope_assim_qsp;
    link->alg = NULL;
    link->state_check = NULL;
    link->CheckConsistency = &CheckConsistency_Nonzero_Model252;
    link->RKSolver = &ExplicitRKSolver;
}


int ReadInitData_Assim_254_qsp(VEC global_params, VEC params, QVSData* qvs, unsigned short int dam, VEC y_0, unsigned int type, unsigned int diff_start, unsigned int no_init_start, void* user, void* external)
{
    //For this type, all initial conditions for variational equation must be set here.
    unsigned int i;
    unsigned int offset = 4;

    y_0.ve[offset++] = 1.0;  //ds_p/ds_p0

    y_0.ve[offset++] = 0.0;  //ds_t/ds_p0

    y_0.ve[offset++] = 0.0;  //ds_s/ds_p0

    y_0.ve[offset++] = 1.0;  //dq/dq_0
    y_0.ve[offset++] = 0.0;  //dq/ds_p0

    for (i = offset; i < y_0.dim; i++)	y_0.ve[i] = 0.0;	//From upstreams

    return 0;
}


//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: A_i,L_i,A_h | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2       3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10     
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void TopLayerHillslope_assim_qsp(double t, VEC y_i, VEC* y_p, unsigned short int num_parents, VEC global_params, double* forcing_values, QVSData* qvs, VEC params, int state, void* user, VEC ans)
{
    unsigned int i, j;

    UpstreamData *updata = (UpstreamData*)user;

    double lambda_1 = global_params.ve[1];
    double k_3 = global_params.ve[4];	//[1/min]
    double h_b = global_params.ve[6];	//[m]
    double S_L = global_params.ve[7];	//[m]
    double A = global_params.ve[8];
    double B = global_params.ve[9];
    double exponent = global_params.ve[10];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params.ve[1];	//[m]
    double A_h = params.ve[2];	//[m^2]
    double invtau = params.ve[3];	//[1/min]
    double k_2 = params.ve[4];	//[1/min]
    double k_i = params.ve[5];	//[1/min]
    double c_1 = params.ve[6];
    double c_2 = params.ve[7];

    double q = y_i.ve[0];		//[m^3/s]
    double s_p = y_i.ve[1];	//[m]
    double s_t = y_i.ve[2];	//[m]
    double s_s = y_i.ve[3];	//[m]

    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    //A few calculations...
    double q_to_lambda_1 = pow(q, lambda_1);
    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);
    double remaining = 1.0 - s_t / S_L;
    double pow_term = (remaining > 0.0) ? pow(remaining, exponent) : 0.0;
    double pow_term_m1 = (remaining > 1e-12) ? pow_term / remaining : pow(1e-12, exponent - 1.0);
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

    //Discharge
    double inflow = 0.0;
    ans.ve[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i < num_parents; i++)
        inflow += y_p[i].ve[0];
    ans.ve[0] = invtau * q_to_lambda_1 * (inflow + ans.ve[0]);

    //Hillslope
    ans.ve[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans.ve[2] = q_pt - q_ts - e_t;
    ans.ve[3] = q_ts - q_sl - e_s;


    //Init for variational equations
    unsigned int offset = 4, dim = ans.dim, problem_dim = 4;
    for (i = offset; i < dim; i++)	ans.ve[i] = 0.0;

    //Compute partial derivatives (local variables)
    double dfq_dq = lambda_1 * invtau * q_to_lambda_1_m1 * (-q + c_2*(k_2*s_p + k_3*s_s) + inflow) - invtau * q_to_lambda_1;
    double dfq_dsp = invtau*q_to_lambda_1*c_2*k_2;
    double dfq_dss = invtau*q_to_lambda_1*c_2*k_3;
    double dfsp_dsp = -k_2 - k_t;
    double dfsp_dst = k_2 / S_L*B*exponent*pow_term_m1*s_p;
    double dfst_dsp = k_t;
    double dfst_dst = -dfsp_dst - k_i;
    double dfss_dst = k_i;
    double dfss_dss = -k_3;

    //Compute partial derivatives (upstreams variables)
    double dfq_dupq = invtau*q_to_lambda_1;

    //Hillslope variational eqs
    ans.ve[offset] = dfsp_dsp * y_i.ve[offset] + dfsp_dst * y_i.ve[offset + 1];	//s_p, s_p
    ans.ve[offset + 1] = dfst_dsp * y_i.ve[offset] + dfst_dst * y_i.ve[offset + 1];	//s_t, s_p
    ans.ve[offset + 2] = dfss_dst * y_i.ve[offset + 1] + dfss_dss * y_i.ve[offset + 2];	//s_s, s_p

    //Discharge variational eqs from this link
    ans.ve[offset + 3] = dfq_dq * y_i.ve[offset + 3]; //q, q
    ans.ve[offset + 4] = dfq_dq * y_i.ve[offset + 4] + dfq_dsp * y_i.ve[offset + 0] + dfq_dss * y_i.ve[offset + 2]; //q, s_p

    //Discharge variational eqs from parent links
    unsigned int current_idx = offset + 5, parent_idx;
    for (i = 0; i < num_parents; i++)
    {
        parent_idx = offset + 3;
        unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
        for (j = 0; j < num_upstreams; j++)
        {
            ans.ve[current_idx] = dfq_dupq * y_p[i].ve[parent_idx] + dfq_dq * y_i.ve[current_idx]; //q, upq
            ans.ve[current_idx + 1] = dfq_dupq * y_p[i].ve[parent_idx + 1] + dfq_dq * y_i.ve[current_idx + 1]; //q, ups_p
            current_idx += 2;
            parent_idx += 2;
        }
    }
}

//I think these values need to be set only for assigned links, not getting.
//fit_states[i] holds the index in each state vector of the ith sensitivity at this link.
//fit_to_universal[i] holds universal index of the ith sensitivity at this link.
//These only store those sensitivites used for the fitting.
void Setup_Fitting_Data_Model254_qsp(AsynchSolver* asynch, unsigned int* obs_locs, unsigned int num_obs)
{
    unsigned int i, j, k, my_N = asynch->my_N, *my_sys = asynch->my_sys, *assignments = asynch->assignments;
    Link *sys = asynch->sys, *current;
    UpstreamData *updata;

    //Number of states to fit
    unsigned int counter;
    unsigned int problem_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int assim_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int num_change_states = 2;	//For q and s_p

    //Find links upstreams from gauges
    bool *is_above_gauges;
    unsigned int *above_gauges;
    unsigned int num_above = GaugeDownstream(asynch, obs_locs, num_obs, &above_gauges, &is_above_gauges);

    //Calculate the number of states needed for the fitting
    //unsigned int allstates_needed = num_above * 2;	//For q and s_p

    for (i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == my_rank)
        {
            current = &sys[obs_locs[i]];
            updata = (UpstreamData*)current->user;

            updata->num_fit_states = 1 + updata->num_upstreams;
            //for(j=0;j<current->num_parents;j++)
            //	updata->num_fit_states += updata->num_upstreams[j];
            updata->num_fit_states *= num_change_states;

            updata->fit_states = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));
            updata->fit_to_universal = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));

            //For q, s_p
            updata->fit_states[0] = problem_dim + 3;	//q, q here
            updata->fit_to_universal[0] = current->location * assim_dim;
            updata->fit_states[1] = problem_dim + 4;	//q, s_p here
            updata->fit_to_universal[1] = current->location * assim_dim + 1;
            counter = 2;
            for (j = 0; j < current->num_parents; j++)
            {
                unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
                Link **upstreams = ((UpstreamData *)updata->upstreams[i]->user)->upstreams;

                for (k = 0; k < num_upstreams; k++)
                {
                    updata->fit_states[counter] = problem_dim + 3 + counter;
                    updata->fit_to_universal[counter] = upstreams[k]->location * assim_dim;
                    updata->fit_states[counter + 1] = problem_dim + 4 + counter;
                    updata->fit_to_universal[counter + 1] = upstreams[k]->location * assim_dim + 1;
                    counter += 2;
                }
            }
            /*
            printf("ID = %u\n",current->ID);
            for(k=0;k<current->num_parents;k++)
                printf("upstreams = %u\n",updata->num_upstreams[k]);
            for(k=0;k<updata->num_fit_states;k++)
            {
            printf("%u ",updata->fit_to_universal[k]);
            }
            printf("\n++++++\n");
            */
        }

    }


    free(is_above_gauges);
    free(above_gauges);
    //return allstates_needed;
}



//Data assimilation model (Model 254, q and s_t) ************************************************************************************


//For modifying q and s_t
void InitRoutines_Assim_254_qst(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    UpstreamData* updata = (UpstreamData*)(link->user);
    unsigned int i, problem_dim = 4;	//Number of model eqs

    link->dim = problem_dim + 5	    //Model eqs + variational eqs from this link
        + updata->num_upstreams * 2;  //Variational eqs from upstreams
    //for(i=0;i<link->num_parents;i++)
    //	link->dim += updata->num_upstreams[i] * 2;	//Variational eqs from upstreams
    link->no_ini_start = 4;
    link->diff_start = 0;

    link->num_dense = link->dim - 3;	//Only q, variational eqs	!!!! Do all of the variational eqs really need to be passed down? !!!!
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;
    for (i = 1; i < link->num_dense; i++)	link->dense_indices[i] = i + 3;

    link->f = &TopLayerHillslope_assim_qst;
    link->alg = NULL;
    link->state_check = NULL;
    link->CheckConsistency = &CheckConsistency_Nonzero_Model252_st;
    link->RKSolver = &ExplicitRKSolver;
}


int ReadInitData_Assim_254_qst(VEC global_params, VEC params, QVSData* qvs, unsigned short int dam, VEC y_0, unsigned int type, unsigned int diff_start, unsigned int no_init_start, void* user, void* external)
{
    //For this type, all initial conditions for variational equation must be set here.
    unsigned int i;
    unsigned int offset = 4;

    y_0.ve[offset++] = 0.0;  //ds_p/ds_t0

    y_0.ve[offset++] = 1.0;  //ds_t/ds_t0

    y_0.ve[offset++] = 0.0;  //ds_s/ds_t0

    y_0.ve[offset++] = 1.0;  //dq/dq_0
    y_0.ve[offset++] = 0.0;  //dq/ds_t0

    for (i = offset; i < y_0.dim; i++)	y_0.ve[i] = 0.0;	//From upstreams

    return 0;
}


//Function for river system with data assimilation. Uses model 252/254.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: A_i,L_i,A_h | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2       3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10     
//y_i[0] = q, y_i[1] = s_p, y_i[2] = s_t, y_i[3] = s_s followed by N entries for the variational equation
//!!!! Note: this actually works out to be the same as the function for qsp, I think... !!!!
void TopLayerHillslope_assim_qst(double t, VEC y_i, VEC* y_p, unsigned short int num_parents, VEC global_params, double* forcing_values, QVSData* qvs, VEC params, int state, void* user, VEC ans)
{
    unsigned int i, j;

    UpstreamData *updata = (UpstreamData*)user;

    double lambda_1 = global_params.ve[1];
    double k_3 = global_params.ve[4];	//[1/min]
    double h_b = global_params.ve[6];	//[m]
    double S_L = global_params.ve[7];	//[m]
    double A = global_params.ve[8];
    double B = global_params.ve[9];
    double exponent = global_params.ve[10];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params.ve[1];	//[m]
    double A_h = params.ve[2];	//[m^2]
    double invtau = params.ve[3];	//[1/min]
    double k_2 = params.ve[4];	//[1/min]
    double k_i = params.ve[5];	//[1/min]
    double c_1 = params.ve[6];
    double c_2 = params.ve[7];

    double q = y_i.ve[0];		//[m^3/s]
    double s_p = y_i.ve[1];	//[m]
    double s_t = y_i.ve[2];	//[m]
    double s_s = y_i.ve[3];	//[m]

    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    //A few calculations...
    double q_to_lambda_1 = pow(q, lambda_1);
    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);
    double remaining = 1.0 - s_t / S_L;
    double pow_term = (remaining > 0.0) ? pow(remaining, exponent) : 0.0;
    double pow_term_m1 = (remaining > 1e-12) ? pow_term / remaining : pow(1e-12, exponent - 1.0);
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

    //Discharge
    double inflow = 0.0;
    ans.ve[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i < num_parents; i++)
        inflow += y_p[i].ve[0];
    ans.ve[0] = invtau * q_to_lambda_1 * (inflow + ans.ve[0]);

    //Hillslope
    ans.ve[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans.ve[2] = q_pt - q_ts - e_t;
    ans.ve[3] = q_ts - q_sl - e_s;


    //Init for variational equations
    unsigned int offset = 4, dim = ans.dim, problem_dim = 4;
    for (i = offset; i < dim; i++)	ans.ve[i] = 0.0;

    //Compute partial derivatives (local variables)
    double dfq_dq = lambda_1 * invtau * q_to_lambda_1_m1 * (-q + c_2*(k_2*s_p + k_3*s_s) + inflow) - invtau * q_to_lambda_1;
    double dfq_dsp = invtau*q_to_lambda_1*c_2*k_2;
    double dfq_dss = invtau*q_to_lambda_1*c_2*k_3;
    double dfsp_dsp = -k_2 - k_t;
    double dfsp_dst = k_2 / S_L*B*exponent*pow_term_m1*s_p;
    double dfst_dsp = k_t;
    double dfst_dst = -dfsp_dst - k_i;
    double dfss_dst = k_i;
    double dfss_dss = -k_3;

    //Compute partial derivatives (upstreams variables)
    double dfq_dupq = invtau*q_to_lambda_1;

    //Hillslope variational eqs
    ans.ve[offset] = dfsp_dsp * y_i.ve[offset] + dfsp_dst * y_i.ve[offset + 1];	//s_p, s_t
    ans.ve[offset + 1] = dfst_dsp * y_i.ve[offset] + dfst_dst * y_i.ve[offset + 1];	//s_t, s_t
    ans.ve[offset + 2] = dfss_dst * y_i.ve[offset + 1] + dfss_dss * y_i.ve[offset + 2];	//s_s, s_t

    //Discharge variational eqs from this link
    ans.ve[offset + 3] = dfq_dq * y_i.ve[offset + 3]; //q, q
    ans.ve[offset + 4] = dfq_dq * y_i.ve[offset + 4] + dfq_dsp * y_i.ve[offset + 0] + dfq_dss * y_i.ve[offset + 2]; //q, s_t

    //Discharge variational eqs from parent links
    unsigned int current_idx = offset + 5, parent_idx;
    for (i = 0; i < num_parents; i++)
    {
        parent_idx = offset + 3;
        unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
        for (j = 0; j < num_upstreams; j++)
        {
            ans.ve[current_idx] = dfq_dupq * y_p[i].ve[parent_idx] + dfq_dq * y_i.ve[current_idx]; //q, upq
            ans.ve[current_idx + 1] = dfq_dupq * y_p[i].ve[parent_idx + 1] + dfq_dq * y_i.ve[current_idx + 1]; //q, ups_t
            current_idx += 2;
            parent_idx += 2;
        }
    }
}

//I think these values need to be set only for assigned links, not getting.
//fit_states[i] holds the index in each state vector of the ith sensitivity at this link.
//fit_to_universal[i] holds universal index of the ith sensitivity at this link.
//These only store those sensitivites used for the fitting.
void Setup_Fitting_Data_Model254_qst(AsynchSolver* asynch, unsigned int* obs_locs, unsigned int num_obs)
{
    unsigned int i, j, k, my_N = asynch->my_N, *my_sys = asynch->my_sys, *assignments = asynch->assignments;
    Link *sys = asynch->sys, *current;
    UpstreamData *updata;

    //Number of states to fit
    unsigned int counter;
    unsigned int problem_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int assim_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int num_change_states = 2;	//For q and s_t

    //Find links upstreams from gauges
    bool *is_above_gauges;
    unsigned int *above_gauges;
    unsigned int num_above = GaugeDownstream(asynch, obs_locs, num_obs, &above_gauges, &is_above_gauges);

    for (i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == my_rank)
        {
            current = &sys[obs_locs[i]];
            updata = (UpstreamData*)current->user;

            updata->num_fit_states = 1 + updata->num_upstreams;
            //for(j=0;j<current->num_parents;j++)
            //	updata->num_fit_states += updata->num_upstreams[j];
            updata->num_fit_states *= num_change_states;

            updata->fit_states = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));
            updata->fit_to_universal = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));

            //For q, s_t
            updata->fit_states[0] = problem_dim + 3;	//q, q here
            updata->fit_to_universal[0] = current->location * assim_dim;
            updata->fit_states[1] = problem_dim + 4;	//q, s_t here
            updata->fit_to_universal[1] = current->location * assim_dim + 2;
            counter = 2;
            for (j = 0; j < current->num_parents; j++)
            {
                unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
                Link **upstreams = ((UpstreamData *)updata->upstreams[i]->user)->upstreams;

                for (k = 0; k < num_upstreams; k++)
                {
                    updata->fit_states[counter] = problem_dim + 3 + counter;
                    updata->fit_to_universal[counter] = upstreams[k]->location * assim_dim;
                    updata->fit_states[counter + 1] = problem_dim + 4 + counter;
                    updata->fit_to_universal[counter + 1] = upstreams[k]->location * assim_dim + 2;
                    counter += 2;
                }
            }
            /*
            printf("ID = %u\n",current->ID);
            //for(k=0;k<current->num_parents;k++)
            //	printf("upstreams = %u\n",updata->num_upstreams[k]);
            for(k=0;k<updata->num_fit_states;k++)
            {
            printf("%u %u\n",updata->fit_states[k],updata->fit_to_universal[k]);
            }
            printf("\n++++++\n");
            */
        }

    }


    free(is_above_gauges);
    free(above_gauges);
    //return allstates_needed;
}

void CheckConsistency_Nonzero_Model252_st(VEC y, VEC params, VEC global_params)
{
    unsigned int i, problem_dim = 4;

    if (y.ve[0] < 1e-14)	y.ve[0] = 1e-14;
    if (y.ve[1] > global_params.ve[7])		y.ve[1] = global_params.ve[7];
    for (i = 1; i < problem_dim; i++)
        if (y.ve[i] < 0.0)	y.ve[i] = 0.0;
}

