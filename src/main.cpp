#include <stdio.h>
#include <mpi.h>
#include <sstream>
#include <vector>
#include "sparse_matrix.h"
#include "mpimatrix.h"
#include "generator.h"
#include "dense_matrix.h"
#include "cg.h"

using namespace std;

void printHelp(int rank)
{
    if (rank != 0) return;
    printf("This program is using MPI to do operations on sparse matrices\n");
    printf("Usage:\n");
    printf("./program <matrix_file_path> <matrix_file_path> <operation>\n");
}

void printLU(int rank, sparse_matrix &L, sparse_matrix &U, const char* name)
{
    if (rank != 0) return;
    printf("%s L:\n", name);
    L.printSparse();
    printf("%s U:\n", name);
    U.printSparse();
}

void printM(int rank, sparse_matrix &M, const char* name)
{
    if (rank != 0) return;
    printf("%s:\n", name);
    M.printSparse();
}

void genRandomM(int rank, int size)
{
    if (rank != 0) return;
    auto randM = Generator(rank, size).GenerateRandomMatrix(10, 10, 25, column_wise);
    printf("Random sparse matrix: \n");
    randM.printSparse();
    auto densM = dense_matrix(randM);
    printf("dense matrix: \n");
    densM.printDense();
}

int main (int argc, char *argv[])
{
    int rank, size;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    MpiMatrixHelper mpi_helper(rank, size);

//        genRandomM(rank, size);
    printHelp(rank);

    auto m1 = mpi_helper.load("/tmp/mat", sparse);
    auto m2 = mpi_helper.load("/tmp/mat", sparse);

    sparse_matrix mult_result = mpi_helper.mul(m1, m2);
    sparse_matrix add_result = mpi_helper.add(m1, m2);

    sparse_matrix L, U, CL, CU;
    mpi_helper.LU(m1, CL, CU);
    mpi_helper.ILU(m1, L, U);

    // CG non-mpi test
    m1.transpose();
    sparse_vector x(4, column_wise);
    x[0] = x[1] = x[2] = x[3] = 0;
    sparse_vector b(4, column_wise);
    b[0] = 1; b[1] = 2; b[2] = 3; b[3] = 4;
    cg(m1, x, b, L, U);
    printf("RESULT X =\n");
    x.print();

    printM(rank, m1, "Matrix 1");
    printM(rank, m2, "Matrix 2");
    printM(rank, mult_result, "Multiplication");
    printM(rank, add_result, "Sum");
    printLU(rank, L, U, "ILU ");
    printLU(rank, CL, CU, "LU ");

    MPI_Finalize();
    return 0;
}
