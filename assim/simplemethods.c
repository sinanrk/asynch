#include "simplemethods.h"

//FAREN(N,X,Y,F,RPAR,IPAR)
/*
void full_system0(unsigned int* bigdim,double* t,double* y,double* result,double* rpar,unsigned long int* ipar)
{
	unsigned int i,j;

	//Unpack the parameters
	Link** sys = (Link**) ipar[0];
	unsigned int N = ipar[1];
	UnivVars* GlobalVars = (UnivVars*) ipar[2];
	double lambda_1 = GlobalVars->global_params->ve[1];

	for(i=0;i<N;i++)
	{
		result[i] = -y[i];
		for(j=0;j<sys[i]->numparents;j++)
			result[i] += y[sys[i]->parents[j]->location];
		result[i] = sys[i]->params->ve[12] * pow(y[i],lambda_1) * result[i];
	}
}
*/

void full_system0(unsigned int* bigdim,double* t,double* y,double* result,double* rpar,unsigned long int* ipar)
{
	unsigned int i,j;

	//Unpack the parameters
	Link** sys = (Link**) ipar[0];
	double invtau;
	UnivVars* GlobalVars = (UnivVars*) ipar[2];
	double lambda_1 = GlobalVars->global_params->ve[1];
/*
int dummy = 1;
while(dummy)
{
	sleep(10);
}
*/
	for(i=0;i<*bigdim;i++)
	{
		invtau = sys[i]->params->ve[12];

		//Discharge
		result[i] = -y[i];
		for(j=0;j<sys[i]->numparents;j++)
			result[i] += y[sys[i]->parents[j]->location];
		result[i] = invtau * pow(y[i],lambda_1) * result[i];
	}
}


//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,c_1,c_2,c_3
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13  14  15 
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
void full_system1(unsigned int* bigdim,double* t,double* y,double* result,double* rpar,unsigned long int* ipar)
{
	unsigned int i,j;

	//Unpack the parameters
	Link** sys = (Link**) ipar[0];
	unsigned int N = ipar[1];
	unsigned int ip1,curr_idx;
	double s53,rain_value,invtau,c_1,c_2,c_3;
	UnivVars* GlobalVars = (UnivVars*) ipar[2];
	double lambda_1 = GlobalVars->global_params->ve[1];

	for(i=0;i<*bigdim;i+=2)
	{
		curr_idx = i/2;
		ip1 = i+1;
		s53 = pow(y[ip1],5.0/3.0);
		invtau = sys[curr_idx]->params->ve[12];
		c_1 = sys[curr_idx]->params->ve[13];
		c_2 = sys[curr_idx]->params->ve[14];
		c_3 = sys[curr_idx]->params->ve[15];

		//Find rainfall	!!!! Can this be done outside the loop? !!!!
		for(j=1;j<sys[curr_idx]->rain->n_times;j++)
			if(sys[curr_idx]->rain->rainfall[j][0] >= *t)	break;
		rain_value = sys[curr_idx]->rain->rainfall[j-1][1];

		//Discharge
		result[i] = -y[i] + c_1 * s53;
		for(j=0;j<sys[curr_idx]->numparents;j++)
			result[i] += y[sys[curr_idx]->parents[j]->location * 2];
		result[i] = invtau * pow(y[i],lambda_1) * result[i];

		//Hillslope
		result[ip1] = c_2 * rain_value - c_3 * s53;
	}
}

void full_system2(unsigned int* bigdim,double* t,double* y,double* result,double* rpar,unsigned long int* ipar)
{
	unsigned int i,j;

	//Unpack the parameters
	Link** sys = (Link**) ipar[0];
	unsigned int N = ipar[1];
	unsigned int ip1,curr_idx;
	double s53;
	UnivVars* GlobalVars = (UnivVars*) ipar[2];
	double lambda_1 = GlobalVars->global_params->ve[1];

	for(i=0;i<*bigdim;i+=2)
	{
		curr_idx = i/2;
		ip1 = i+1;
		s53 = pow(y[ip1],5.0/3.0);

		result[i] = -y[i] + sys[curr_idx]->params->ve[14] * s53;
		for(j=0;j<sys[curr_idx]->numparents;j++)
			result[i] += y[sys[curr_idx]->parents[j]->location];
		result[i] = sys[curr_idx]->params->ve[12] * pow(y[i],lambda_1) * result[i];

		result[ip1] = -sys[curr_idx]->params->ve[17] * s53;
	}
}

/*
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,RC,u_0
//The numbering is:        0      1        2     3   4   5  6
void full_system200(unsigned int* bigdim,double* t,double* y,double* result,double* rpar,unsigned long int* ipar)
{
	unsigned int i,j;

	//Unpack the parameters
	Link** sys = (Link**) ipar[0];
	Link *current,*child;
	unsigned int N = ipar[1];
	UnivVars* GlobalVars = (UnivVars*) ipar[2];
	//double lambda_1 = GlobalVars->global_params->ve[1];
	double u_0 = GlobalVars->global_params->ve[6];
	double v_0 = GlobalVars->global_params->ve[0];
	double l,qd,S,parent_sum;
	unsigned int ip1,p_loc;
	double loss_coef = u_0/v_0 + 1.0;

	//y[i/2] contains the downstream discharge for link i
	//y[i/2+1] contains the storage for link i
	for(i=0;i<*bigdim;i+=2)
	{
		current = sys[i/2];
		ip1 = i+1;
		l = current->params->ve[0];
		qd = y[i];
		S = y[ip1];

		//Loss
		if(current->numparents > 0)
			result[ip1] = -loss_coef * qd;
		else
			result[ip1] = -qd;

		//Gain
		for(j=0;j<current->numparents;j++)
		{
			p_loc = current->parents[j]->location;
			result[ip1] += y[2*p_loc];
		}
		if(child != NULL)
		{
			parent_sum = 0.0;
			for(j=0;j<child->numparents;j++)
				parent_sum += y[2*child->parents[j]->location];
			result[ip1] += qd / parent_sum * u_0/v_0 * y[2*child->location];
		}


		//Convert to m^3 / min
		result[ip1] *= 60.0;

		//Calculate derivative of qd
		result[i] = v_0/l * result[ip1];
	}
}
*/


//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,RC,u_0
//The numbering is:        0      1        2     3   4   5  6
void full_system200(unsigned int* bigdim,double* t,double* y,double* result,double* rpar,unsigned long int* ipar)
{
	unsigned int i,j;

	//Unpack the parameters
	Link** sys = (Link**) ipar[0];
	Link *current,*child;
	unsigned int N = ipar[1];
	UnivVars* GlobalVars = (UnivVars*) ipar[2];
	//double lambda_1 = GlobalVars->global_params->ve[1];
	double u_0 = GlobalVars->global_params->ve[6];
	double v_0 = GlobalVars->global_params->ve[0];
	double l,qd,S,parent_sum;
	unsigned int ip1,p_loc;

	//y[i/2] contains the downstream discharge for link i
	//y[i/2+1] contains the storage for link i
	for(i=0;i<*bigdim;i+=2)
	{
		current = sys[i/2];
		ip1 = i+1;
		l = current->params->ve[0];
		qd = y[i];
		S = y[ip1];

		//Child
		result[ip1] = -qd;
		child = current->c;
		if(child != NULL)
		{
			//result[ip1] += current->params->ve[2] / (child->params->ve[2] - child->params->ve[1]) * u_0/child->params->ve[0] * y[2*child->location+1];

			parent_sum = 0.0;
			for(j=0;j<child->numparents;j++)
				parent_sum += y[2*child->parents[j]->location];
			if(parent_sum > 1e-12)
				result[ip1] += qd / parent_sum * u_0/child->params->ve[0] * y[2*child->location+1];
		}
		
		//Parents
		if(current->numparents > 0)
		{
			result[ip1] += -u_0/l * S;
			for(j=0;j<current->numparents;j++)
			{
				p_loc = current->parents[j]->location;
				result[ip1] += y[2*p_loc];
			}
		}

		result[ip1] *= 60.0;
		result[i] = (u_0+v_0)/l * result[ip1];
	}
}

