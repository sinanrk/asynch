#include <stdlib.h>
#include "asynch_interface.h"
#include "assim_ls_methods.h"

int* Partition_METIS_ByEqs(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting);

void Setup_Errors(asynchsolver* asynch,unsigned int problem_dim);
unsigned int BuildStateShift(asynchsolver* asynch,unsigned int allstates,unsigned int* data_locs,unsigned int numdata,unsigned int** vareq_shift,unsigned int** inv_vareq_shift);

//Assim model 15
void SetParamSizes_Assim(UnivVars* GlobalVars,void* external);
void ConvertParams_Assim(VEC params,unsigned int type,void* external);
void InitRoutines_Assim(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
void InitRoutines_Model(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
void Precalculations_Assim(Link* link_i,VEC global_params,VEC params,unsigned int disk_params,unsigned int params_size,unsigned short int dam,unsigned int type,void* external);
int ReadInitData_Assim(VEC global_params,VEC params,QVSData* qvs,unsigned short int dam,VEC y_0,unsigned int type,unsigned int diff_start,unsigned int no_init_start,void* user,void* external);
void assim_river_rainfall_adjusted_custom(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);

//Assim model 254
void SetParamSizes_Assim_254(UnivVars* GlobalVars,void* external);
void ConvertParams_Assim_254(VEC params,unsigned int type,void* external);
void InitRoutines_Assim_254(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
void InitRoutines_Model_254(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
void Precalculations_Assim_254(Link* link_i,VEC global_params,VEC params,unsigned int disk_params,unsigned int params_size,unsigned short int dam,unsigned int type,void* external);
int ReadInitData_Assim_254(VEC global_params,VEC params,QVSData* qvs,unsigned short int dam,VEC y_0,unsigned int type,unsigned int diff_start,unsigned int no_init_start,void* user,void* external);
void CheckConsistency_Nonzero_Model254(VEC y,VEC params,VEC global_params);
void TopLayerHillslope_extras_assim(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void Setup_Fitting_Data_Model254(asynchsolver* asynch,unsigned int* data_locs,unsigned int numdata);

//Assim model 254, q
void InitRoutines_Assim_254_q(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
void InitRoutines_Model_252(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
int ReadInitData_Assim_254_q(VEC global_params,VEC params,QVSData* qvs,unsigned short int dam,VEC y_0,unsigned int type,unsigned int diff_start,unsigned int no_init_start,void* user,void* external);
void CheckConsistency_Nonzero_Model252(VEC y,VEC params,VEC global_params);
void TopLayerHillslope_assim_q(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void Setup_Fitting_Data_Model254_q(asynchsolver* asynch,unsigned int* data_locs,unsigned int numdata);

//Assim model 254, q and s_p
void InitRoutines_Assim_254_qsp(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
int ReadInitData_Assim_254_qsp(VEC global_params,VEC params,QVSData* qvs,unsigned short int dam,VEC y_0,unsigned int type,unsigned int diff_start,unsigned int no_init_start,void* user,void* external);
void TopLayerHillslope_assim_qsp(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void Setup_Fitting_Data_Model254_qsp(asynchsolver* asynch,unsigned int* data_locs,unsigned int numdata);

//Assim model 254, q and s_t
void InitRoutines_Assim_254_qst(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
int ReadInitData_Assim_254_qst(VEC global_params,VEC params,QVSData* qvs,unsigned short int dam,VEC y_0,unsigned int type,unsigned int diff_start,unsigned int no_init_start,void* user,void* external);
void TopLayerHillslope_assim_qst(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void Setup_Fitting_Data_Model254_qst(asynchsolver* asynch,unsigned int* data_locs,unsigned int numdata);
void CheckConsistency_Nonzero_Model252_st(VEC y,VEC params,VEC global_params);

