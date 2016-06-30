#include "structs.h"

void full_system0(unsigned int* bigdim,double* t,double* y,double* result,double* rpar,unsigned long int* ipar);
void full_system1(unsigned int* bigdim,double* t,double* y,double* result,double* rpar,unsigned long int* ipar);
void full_system2(unsigned int* bigdim,double* t,double* y,double* result,double* rpar,unsigned long int* ipar);
void full_system200(unsigned int* bigdim,double* t,double* y,double* result,double* rpar,unsigned long int* ipar);

//Fortran functions called from C code
void dopri5_(unsigned int* bigdim,void (*full_system) (unsigned int*,double*,double*,double*,double*,unsigned long int*),double* t_0,double* y_0,double* t_f,double* rel,double* abs,int* itol,void solout(void*,void*,void*,void*,void*,void*,void*,void*,void*,void*,void*),int* iout,double* work,unsigned int* lwork,int* iwork,unsigned int* liwork,double* rpar,unsigned long int* ipar,int* idid);
void openoutputfile_();
void closeoutputfile_();

