//
// Created by ghash on 04.05.15.
//

#ifndef MPI_MATRICES_CG_H
#define MPI_MATRICES_CG_H

#include "sparse_matrix.h"

sparse_vector solveILU(const sparse_matrix &L, const sparse_matrix &U, const sparse_vector &b);
int cg(const sparse_matrix &A, sparse_vector &x, sparse_vector &b, sparse_matrix &L, sparse_matrix &U);

#endif //MPI_MATRICES_CG_H
