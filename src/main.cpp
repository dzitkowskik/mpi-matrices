#include <stdio.h>
#include <mpi.h>
#include <sstream>
#include <vector>
#include "sparse_matrix.h"
#include "mpimatrix.h"

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

    if (rank == 0) printHelp();

    MpiMatrix m1(rank, size);
    MpiMatrix m2(rank, size);

    m1.loadFromFile("/tmp/matrix1", column_wise);
    m2.loadFromFile("/tmp/matrix2", row_wise);

    MpiMatrix mult_result = m1*m2;
    MpiMatrix add_result = m1+m2;

    MpiMatrix L, U;
    m1.LU(L, U);

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
    }
    MPI_Finalize();
    return 0;
}
