#ifndef MPI_MATRIX_H
#define MPI_MATRIX_H

#include <mpi.h>
#include "sparse_matrix.h"

enum MatrixType { sparse, dense, MatrixType_count };

class MpiMatrixHelper
{
public:
	int rank;
	int processors_cnt;
	MPI_Datatype sparse_elem_type;

public:
	MpiMatrixHelper();
	MpiMatrixHelper(int rank, int proc_cnt);
	MpiMatrixHelper(const MpiMatrixHelper &m);
	~MpiMatrixHelper();

public:
	sparse_matrix load(const char* path, MatrixType type, direction dir = column_wise);

	sparse_matrix add(const sparse_matrix &a, const sparse_matrix &b);
	sparse_matrix sub(const sparse_matrix &a, const sparse_matrix &b);
	sparse_matrix mul(const sparse_matrix &a, const sparse_matrix &b);
//	sparse_matrix div(const sparse_matrix &a, const sparse_matrix &b);
	void addto(sparse_matrix &to, const sparse_matrix &what);
	void subto(sparse_matrix &to, const sparse_matrix &what);
//	void multo(sparse_matrix &to, const sparse_matrix &what);
//	void divto(sparse_matrix &to, const sparse_matrix &what);

	sparse_vector mul(const sparse_matrix &A, const sparse_vector &x);

	void LU(const sparse_matrix &A, sparse_matrix &L, sparse_matrix &U);
	void ILU(const sparse_matrix &A, sparse_matrix &L, sparse_matrix &U);
	sparse_vector CG(const sparse_matrix &A, const sparse_vector &b);

private:
	void init();
	void createSparseElemDatatype();
	void sendMatrix(int node, sparse_matrix matrix);
	sparse_matrix receiveMatrix(int node, direction dir);
	void sendVector(int node, sparse_vector vector);
	sparse_vector receiveVector(int node);
};

#endif
