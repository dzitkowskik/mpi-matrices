//
// Created by ghash on 29.03.15.
//

#ifndef MPI_MATRICES_SPARSE_ELEM_H
#define MPI_MATRICES_SPARSE_ELEM_H

struct sparse_matrix_elem
{
	int col;
	int row;
	double value;
};

#endif //MPI_MATRICES_SPARSE_ELEM_H
