#include <stdio.h>
#include <mpi.h>
#include <sstream>
#include <vector>
#include "sparse_matrix.h"
#include "mpimatrix.h"
#include "generator.h"
#include "dense_matrix.h"

using namespace std;

void printHelp()
{
    printf("This program is using MPI to do operations on sparse matrices\n");
    printf("Usage:\n");
    printf("./program <matrix_file_path> <matrix_file_path> <operation>\n");
}

int main (int argc, char *argv[])
{
    int rank, size;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

//    auto randM = Generator(rank, size).GenerateRandomMatrix(10, 10, 25, column_wise);
//    printf("Random sparse matrix: \n");
//    randM.print();
//    auto densM = dense_matrix(randM.matrix);
//    printf("dense matrix: \n");
//    densM.printDense();

    MpiMatrix m1, m2;
    if (rank == 0)
        printHelp();

    m1 = MpiMatrix::load("/tmp/mat", rank, size, sparse);
    m2 = MpiMatrix::load("/tmp/mat", rank, size, sparse);

    MpiMatrix mult_result = m1*m2;
    MpiMatrix add_result = m1+m2;

    MpiMatrix CL, CU;
    m1.LU(CL, CU);

    MpiMatrix L, U;
    m1.ILU(L, U);

    if (rank == 0)
    {
        printf("Matrix 1: \n");
        m1.print();
        printf("Matrix 2: \n");
        m2.print();

        printf("Multiplication:\n");
        mult_result.print();
        printf("Sum:\n");
        add_result.print();
        printf("L:\n");
        L.print();
        printf("U:\n");
        U.print();
        printf("CL:\n");
        CL.print();
        printf("CU:\n");
        CU.print();
    }
    MPI_Finalize();
    return 0;
}
