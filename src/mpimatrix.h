#ifndef MPI_MATRIX_H
#define MPI_MATRIX_H


#include <mpi.h>
#include "sparse_matrix.h"

class MpiMatrix
{
public:
	int rank;
	int processors_cnt;
	sparse_matrix matrix;
	MPI_Datatype sparse_elem_type;

public:
	MpiMatrix();
	MpiMatrix(int rank, int proc_cnt);
	MpiMatrix(int rank, int proc_cnt, sparse_matrix sp);
	MpiMatrix(const MpiMatrix &m);
	~MpiMatrix();

public:
	void print();
	void loadFromFile(const char *path, direction dir);
	MpiMatrix operator+(const MpiMatrix &m);
	MpiMatrix operator*(const MpiMatrix &m);
	void LU(MpiMatrix &L, MpiMatrix &U);

private:
	void init();
	void createSparseElemDatatype();
	void sendMatrix(int node, sparse_matrix matrix);
	sparse_matrix receiveMatrix(int node, direction dir);

	void sendVector(int node, sparse_vector vector);
	sparse_vector receiveVector(int node);
};

#endif
