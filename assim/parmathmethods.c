#include "parmathmethods.h"


extern int flaggy;

//Calculates A = V*D^-1*V^T in parallel, where D is diagonal and square.
//V and D are assumed square.
//This assumes D is on all processes, but V is distributed.
//temp (same size as Vlocal) are for temporary storage.
//Vlocal is assumed transposed.
void parVDinvVT_mlt(MAT* Vlocal,int* descV,VEC* D,MAT* Alocal,MAT* temp,int my_col,int npcol)
{
	unsigned int i,j,globalcol;
	int ia,ja,ib,jb,ic,jc;
	unsigned int locrow = Vlocal->m;
	unsigned int loccol = Vlocal->n;
	unsigned int col_blocking = descV[5];
	double dinv;
	char transa = 'N';	//!!!! May need to play with these to get the right thing !!!!
	char transb = 'T';
	double alpha = 1.0;
	double beta = 0.0;
	ia = ja = ib = jb = ic = jc = 1;

	//First compute V*D^-1
	for(i=0;i<locrow;i++)
	{
		globalcol = my_col*col_blocking + (i/col_blocking)*(col_blocking*npcol) + i%col_blocking;
		dinv = 1.0 / D->ve[globalcol];
		for(j=0;j<loccol;j++)
			temp->me[i][j] = Vlocal->me[i][j] * dinv;
	}

	//Now compute temp*V^T
	pdgemm_(&transa,&transb,&(descV[2]),&(descV[3]),&(descV[3]),&alpha,temp->array,&ia,&ja,descV,Vlocal->array,&ib,&jb,descV,&beta,Alocal->array,&ic,&jc,descV);
}


//Calculates A = V*sqrt(D)^-1*V^T in parallel, where D is diagonal and square.
//V and D are assumed square.
//This assumes D is on all processes, but V is distributed.
//temp (same size as Vlocal) are for temporary storage.
//Vlocal is assumed transposed, and Alocal is given transposed.
void parVsqrtDinvVT_mlt(MAT* Vlocal,int* descV,VEC* D,MAT* Alocal,MAT* temp,int my_col,int npcol)
{
	unsigned int i,j,globalcol;
	int ia,ja,ib,jb,ic,jc;
	unsigned int locrow = Vlocal->m;
	unsigned int loccol = Vlocal->n;
	unsigned int col_blocking = descV[5];
	double sqrtdinv;
	char transa = 'N';	//!!!! May need to play with these to get the right thing !!!!
	char transb = 'T';
	double alpha = 1.0;
	double beta = 0.0;
	ia = ja = ib = jb = ic = jc = 1;

	//First compute V*sqrt(D)^-1
	for(i=0;i<locrow;i++)
	{
		globalcol = my_col*col_blocking + (i/col_blocking)*(col_blocking*npcol) + i%col_blocking;
		sqrtdinv = 1.0 / sqrt(D->ve[globalcol]);
		for(j=0;j<loccol;j++)
			temp->me[i][j] = Vlocal->me[i][j] * sqrtdinv;
	}

	//Now compute temp*V^T
	pdgemm_(&transa,&transb,&(descV[2]),&(descV[3]),&(descV[3]),&alpha,temp->array,&ia,&ja,descV,Vlocal->array,&ib,&jb,descV,&beta,Alocal->array,&ic,&jc,descV);
}


//Computes B = A^T*D, with A and B distributed and D diagonal but given on all procs.
//B has the same descriptor as A.
//A and B are assumed transposed.
void parmTdiag_mlt(MAT* Alocal,int* descA,VEC* D,MAT* Blocal,int my_col,int npcol)
{
	unsigned int i,j,globalcol;
	unsigned int locrow = Blocal->m;
	unsigned int loccol = Blocal->n;
	unsigned int col_blocking = descA[5];

	//Compute A*D
	for(j=0;j<loccol;j++)
	{
		globalcol = my_col*col_blocking + (j/col_blocking)*(col_blocking*npcol) + j%col_blocking;
		for(i=0;i<locrow;i++)
			Blocal->me[j][i] = Alocal->me[i][j] * D->ve[globalcol];
	}
}

//Transposes the distributed matrix A and multiplies by alpha.
//The result is stored as B.
void parTranspose(MAT* Alocal,int* descA,double alpha,MAT* Clocal,int* descC)
{
	int ia,ja,ic,jc;
	double beta = 0.0;
	ia = ja = ic = jc = 1;

	pdtran_(&(descC[2]),&(descC[3]),&alpha,Alocal->array,&ia,&ja,descA,&beta,Clocal->array,&ic,&jc,descC);
}

//Calculates C = A*B in parallel.
//Alocal and Blocal are assumed to be transposed.
void parmm_mlt(MAT* Alocal,int* descA,MAT* Blocal,int* descB,MAT* Clocal,int* descC)
{
	int ia,ja,ib,jb,ic,jc;
	char transa = 'N';	//!!!! May need to play with these to get the right thing !!!!
	char transb = 'N';
	double alpha = 1.0;
	double beta = 0.0;
	ia = ja = ib = jb = ic = jc = 1;

	pdgemm_(&transa,&transb,&(descA[2]),&(descB[3]),&(descA[3]),&alpha,Alocal->array,&ia,&ja,descA,Blocal->array,&ib,&jb,descB,&beta,Clocal->array,&ic,&jc,descC);
	//pdgemm_(&transa,&transb,&(descA[3]),&(descB[2]),&(descA[2]),&alpha,Alocal->array,&ia,&ja,descA,Blocal->array,&ib,&jb,descB,&beta,Clocal->array,&ic,&jc,descC);
}

//Calculates w = A*v
void parmv_mlt(MAT* Alocal,int* descA,VEC* vlocal,int* descv,VEC* wlocal,int* descw)
{
	int ia,ja,iv,jv,iw,jw;
	int inc = 1;
	char T = 'N';
	double alpha = 1.0;
	double beta = 0.0;
	ia = ja = iv = jv = iw = jw = 1;

	pdgemv_(&T,&(descA[2]),&(descA[3]),&alpha,Alocal->array,&ia,&ja,descA,vlocal->ve,&iv,&jv,descv,&inc,&beta,wlocal->ve,&iw,&jw,descw,&inc);
}


//Calculates A = A + alpha*I.
void parmaeye_add(MAT* Alocal,int* descA,double alpha,int my_row,int my_col,int nprow,int npcol)
{
	unsigned int i,j,globalrow,globalcol;
	unsigned int loc_rows = Alocal->m;
	unsigned int loc_cols = Alocal->n;
	unsigned int row_blocking = descA[4];
	unsigned int col_blocking = descA[5];

	for(i=0;i<loc_rows;i++)
	{
		for(j=0;j<loc_cols;j++)
		{
			globalcol = my_row*row_blocking + (j/row_blocking)*(row_blocking*nprow) + j%row_blocking;
			globalrow = my_col*col_blocking + (i/col_blocking)*(col_blocking*npcol) + i%col_blocking;
			if(globalrow == globalcol)
				Alocal->me[i][j] += alpha;
		}
	}
}

//Calculates w = v + alpha*u, with w and v distributed, and u global.
void parvglobal_add(VEC* vlocal,int* descv,double alpha,VEC* u,VEC* wlocal,int my_row,int nprow)
{
	unsigned int i,globalrow;
	unsigned int loc_rows = vlocal->dim;
	unsigned int col_blocking = descv[5];

	for(i=0;i<loc_rows;i++)
	{
		globalrow = my_row*col_blocking + (i/col_blocking)*(col_blocking*nprow) + i%col_blocking;
		wlocal->ve[i] = vlocal->ve[i] + alpha * u->ve[globalrow];
	}
}

//Calculates the average across each row of A.
//Note: This assumes Alocal is tranposed, so it's really going
//to compute across the cols.
//temp has same dim as ave.
void parensemble_average(MAT* Alocal,int* descA,VEC* ave,int my_row,int nprow,unsigned int ensemble_size,VEC* temp,int applyflag)
{
	unsigned int i,j,globalcol;
	unsigned int numcols = ave->dim;
	unsigned int row_blocking = descA[4];
	//unsigned int col_blocking = descA[5];
	unsigned int loc_rows = Alocal->m;
	unsigned int loc_cols = Alocal->n;

	//Initialize ave and temp
	for(i=0;i<numcols;i++)
	{
		ave->ve[i] = 0.0;
		temp->ve[i] = 0.0;
	}

	//Calculate each component of ave on this proc
	for(j=0;j<loc_cols;j++)
	{
		globalcol = my_row*row_blocking + (j/row_blocking)*(row_blocking*nprow) + j%row_blocking;
		for(i=0;i<loc_rows;i++)
			temp->ve[globalcol] += Alocal->me[i][j];
	}

	//Compute average
	for(i=0;i<numcols;i++)	temp->ve[i] = temp->ve[i] / ensemble_size;

	//Send average to each proc
	MPI_Allreduce(temp->ve,ave->ve,numcols,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	if(applyflag)
	{
		//Subtract the average from each row
		for(j=0;j<loc_cols;j++)
		{
			globalcol = my_row*row_blocking + (j/row_blocking)*(row_blocking*nprow) + j%row_blocking;
			for(i=0;i<loc_rows;i++)
				Alocal->me[i][j] -= ave->ve[globalcol];
		}
	}
}

//Add v[i] to every entry in col i of A.
//temp and temp2 should have dimension ensemble_size.
void parapply_vector(MAT* Alocal,int* descA,VEC* vlocal,int* descv,int my_row,int nprow,VEC* temp,VEC* temp2)
{
	unsigned int i,j,globalrow;
	unsigned int m = descA[2];
	unsigned int row_blocking = descA[4];
	unsigned int loc_rows = Alocal->m;
	unsigned int loc_cols = Alocal->n;

	//Make sure proc 0 is in the first row and col
	if(descA[6] != 0 && descA[7] != 0)
	{
		printf("Error in parapply_vector: Need the first proc to be 0.\n");
		return;
	}

	//Make v global !!!! This can probably be done better !!!!
	for(i=0;i<m;i++)
	{
		if( (int)((i/row_blocking) % nprow) == my_row )
			temp->ve[i] = vlocal->ve[(i/(nprow*row_blocking))*row_blocking+i%row_blocking];
		else
			temp->ve[i] = 0.0;
	}
	MPI_Allreduce(temp->ve,temp2->ve,m,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	//Add v from each row
	for(j=0;j<loc_cols;j++)
	{
		globalrow = my_row*row_blocking + (j/row_blocking)*(row_blocking*nprow) + j%row_blocking;
		for(i=0;i<loc_rows;i++)
			Alocal->me[i][j] += temp2->ve[globalrow];
	}
}

//Add v[i] to all every entry in col i of A.
void parapply_globalvector(MAT* Alocal,int* descA,VEC* v,int my_row,int nprow)
{
	unsigned int i,j,globalrow;
	unsigned int row_blocking = descA[4];
	unsigned int loc_rows = Alocal->m;
	unsigned int loc_cols = Alocal->n;

	//Make sure proc 0 is in the first row and col
	if(descA[6] != 0 && descA[7] != 0)
	{
		printf("Error in parapply_vector: Need the first proc to be 0.\n");
		return;
	}

	//Add v from each row
	for(j=0;j<loc_cols;j++)
	{
		globalrow = my_row*row_blocking + (j/row_blocking)*(row_blocking*nprow) + j%row_blocking;
		for(i=0;i<loc_rows;i++)
			Alocal->me[i][j] += v->ve[globalrow];
	}
}

//Makes the Ith row of matrix A fully available on all processes.
void parglobalize_row(MAT* Alocal,int* descA,int my_row,int my_col,int nprow,int npcol,unsigned int rowidx,double* storage)
{
	unsigned int i,j,size,curr_proc;
	unsigned int m = descA[2];
	unsigned int row_blocking = descA[4];
	unsigned int col_blocking = descA[5];

	///Find processor row
	unsigned int proc_row = (rowidx/col_blocking)%npcol;
	unsigned int proc_col = 0;
	unsigned int start_proc = proc_col * npcol + proc_row;	//Switched proc_row and proc_col
	curr_proc = start_proc;
	unsigned int locrowidx = rowidx % col_blocking + rowidx/(npcol*col_blocking)*col_blocking;
	unsigned int loccolidx;

	for(j=0;j<m;j+=col_blocking)
	{
		//Broadcast the data
		size = (j+col_blocking < m) ? col_blocking : m - j;
		if((int)curr_proc == my_rank)
		{
			loccolidx = j % row_blocking + j/(nprow*row_blocking)*row_blocking;
			for(i=loccolidx;i<loccolidx+size;i++)
				storage[j+i-loccolidx] = Alocal->me[locrowidx][i];
		}

		MPI_Bcast(&(storage[j]),size,MPI_DOUBLE,curr_proc,MPI_COMM_WORLD);
		proc_col = (proc_col + 1) % nprow;
		curr_proc = proc_col * npcol + proc_row;
	}

	MPI_Barrier(MPI_COMM_WORLD);

}

//Makes storage the Ith row of matrix A.
void parlocallize_row(MAT* Alocal,int* descA,int my_row,int my_col,int nprow,int npcol,unsigned int rowidx,double* storage)
{
	unsigned int i,j,size,curr_proc;
	unsigned int m = descA[2];
	unsigned int col_blocking = descA[5];

	///Find processor row
	unsigned int proc_row = (rowidx/col_blocking)%npcol;
	unsigned int proc_col = 0;
	unsigned int start_proc = proc_col * npcol + proc_row;	//Switched proc_row and proc_col
	curr_proc = start_proc;
	unsigned int locrowidx = rowidx % col_blocking + rowidx/(npcol*col_blocking)*col_blocking;
	unsigned next_place = 0;


	for(j=0;j<m;j+=col_blocking)
	{
		if( (int)curr_proc == my_rank)
		{
			size = (j+col_blocking < m) ? col_blocking : m - j;
			for(i=next_place;i<next_place+size;i++)
			{
				Alocal->me[locrowidx][i] = storage[j+i-next_place];
			}
			next_place = i;
		}

		proc_col = (proc_col + 1) % nprow;
		curr_proc = proc_col * npcol + proc_row;
	}

	MPI_Barrier(MPI_COMM_WORLD);
}


//Prints the distributed matrix A to stdout.
void parPrint_Matrix(MAT* Alocal,int* descA,int my_row,int my_col,int nprow,int npcol)
{
	unsigned int i,j;
	unsigned int row_blocking = descA[4];
	unsigned int col_blocking = descA[5];
	unsigned int m = descA[2];
	unsigned int n = descA[3];
	MAT* temp = m_get(m,n);
	MAT* temp2 = m_get(m,n);

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			if( (int)((i/row_blocking) % nprow) == my_row && (int)((j/col_blocking) % npcol) == my_col )
				temp->me[i][j] = Alocal->me[(j/(npcol*col_blocking))*col_blocking+j%col_blocking][(i/(nprow*row_blocking))*row_blocking+i%row_blocking];
				//temp->me[i][j] = Alocal->me[(i/(nprow*row_blocking))*row_blocking+i%row_blocking][(j/(npcol*col_blocking))*col_blocking+j%col_blocking];
		}
	}

	MPI_Reduce(temp->array,temp2->array,m*n,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	if(my_rank == 0)	Print_Matrix(temp2);
	m_free(temp);
	m_free(temp2);
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(1);
}

//Prints the distributed vector A to stdout.
void parPrint_Vector(VEC* vlocal,int* descv,int my_row,int my_col,int nprow,int npcol)
{
	unsigned int i;
	unsigned int row_blocking = descv[4];
	unsigned int m = descv[2];
	VEC* temp = v_get(m);
	VEC* temp2 = v_get(m);

	for(i=0;i<m;i++)
	{
		if( (int)((i/row_blocking) % nprow) == my_row )
			temp->ve[i] = vlocal->ve[(i/(nprow*row_blocking))*row_blocking+i%row_blocking];
	}

	MPI_Reduce(temp->ve,temp2->ve,m,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	if(my_rank == 0)	Print_Vector(temp2);
	v_free(temp);
	v_free(temp2);
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(1);
}

