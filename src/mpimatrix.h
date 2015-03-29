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

	MpiMatrix()
	{ init(); }

	MpiMatrix(int rank, int proc_cnt)
			: rank(rank), processors_cnt(proc_cnt)
	{ init(); }

	MpiMatrix(int rank, int proc_cnt, sparse_matrix sp)
			: rank(rank), processors_cnt(proc_cnt), matrix(sp)
	{ init(); }

	MpiMatrix(const MpiMatrix &m)
	{
		matrix = m.matrix;
		rank = m.rank;
		processors_cnt = m.processors_cnt;
		sparse_elem_type = m.sparse_elem_type;
	}

	~MpiMatrix()
	{ }

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

	MPI_Datatype sparse_elem_type;
};

#endif
