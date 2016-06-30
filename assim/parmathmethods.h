#ifndef PARMATHMETHODS_H
#define PARMATHMETHODS_H

#include "mathmethods.h"
#include "my_scalapack.h"
#include "mpi.h"
#include <unistd.h>

extern int my_rank;
extern int np;

void parVDinvVT_mlt(MAT* Vlocal,int* descV,VEC* D,MAT* Alocal,MAT* temp,int my_col,int npcol);
void parVsqrtDinvVT_mlt(MAT* Vlocal,int* descV,VEC* D,MAT* Alocal,MAT* temp,int my_col,int npcol);
void parmTdiag_mlt(MAT* Alocal,int* descA,VEC* D,MAT* Blocal,int my_col,int npcol);
void parTranspose(MAT* Alocal,int* descA,double alpha,MAT* Clocal,int* descC);
void parmm_mlt(MAT* Alocal,int* descA,MAT* Blocal,int* descB,MAT* Clocal,int* descC);
void parmv_mlt(MAT* Alocal,int* descA,VEC* vlocal,int* descv,VEC* wlocal,int* descw);
void parmaeye_add(MAT* Alocal,int* descA,double alpha,int my_row,int my_col,int nprow,int npcol);
void parvglobal_add(VEC* vlocal,int* descv,double alpha,VEC* u,VEC* wlocal,int my_col,int npcol);
void parensemble_average(MAT* Alocal,int* descA,VEC* ave,int my_row,int nprow,unsigned int ensemble_size,VEC* temp,int applyflag);
void parapply_vector(MAT* Alocal,int* descA,VEC* vlocal,int* descv,int my_row,int nprow,VEC* temp,VEC* temp2);
void parapply_globalvector(MAT* Alocal,int* descA,VEC* v,int my_row,int nprow);
void parglobalize_row(MAT* Alocal,int* descA,int my_row,int my_col,int nprow,int npcol,unsigned int rowidx,double* storage);
void parlocallize_row(MAT* Alocal,int* descA,int my_row,int my_col,int nprow,int npcol,unsigned int rowidx,double* storage);

void parPrint_Matrix(MAT* Alocal,int* descA,int my_row,int my_col,int nprow,int npcol);
void parPrint_Vector(VEC* vlocal,int* descv,int my_row,int my_col,int nprow,int npcol);

#endif

